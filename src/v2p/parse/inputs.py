"""Parsers that normalise the four HCC1395 evidence types into one schema.

Every parser yields dicts with a common core so that stage 2 can iterate
over a single manifest:

    variant_id      stable, deterministic, human-readable
    variant_class   SNV | INDEL | MNV | RNA_EDITING | FUSION | AS_<TYPE>
    source          logical name of the input file
    genes           semicolon-joined gene symbols / ids
    locus           display string for the genomic location
    confidence      numeric or '' - meaning differs per source (documented)
    payload         source-specific fields kept verbatim for stage 2
"""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Any, Iterator

# --------------------------------------------------------------------------
# 1. Somatic SNV / InDel VCF  (SEQC2 high-confidence sSNV/sIndel call set)
# --------------------------------------------------------------------------

def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "rt", encoding="utf-8", errors="replace")


def classify_small_variant(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "INDEL"


_AF_KEYS = ("AF", "VAF", "AF_ALT")


def _info_af(info: str) -> float | None:
    """Allele frequency from an INFO field, or None if it says nothing.

    `AF` is the VCF specification's key. `VAF` is what several somatic
    callers write instead - the SEQC2 truth set among them - and it
    means the same thing here, so both are read. A multi-allelic value
    ("0.2,0.8") takes the first, matching the row's first ALT.

    None means "not stated", which is not the same as zero: a variant
    whose frequency is unknown must not be silently dropped by a
    frequency filter.
    """
    if not info:
        return None
    for field in info.split(";"):
        key, _, value = field.partition("=")
        if key in _AF_KEYS and value:
            try:
                return float(value.split(",")[0])
            except ValueError:
                return None
    return None


def _phase(fields: list[str], alt_index: int) -> tuple[str, frozenset]:
    """(phase set, haplotypes carrying this ALT) from the first sample.

    Only a pipe-separated GT is phase. `0/1` is a genotype: it says the
    sample is heterozygous, not which allele the variant sits on, so it
    yields nothing here. `PS` names the phase block, and two variants
    share a haplotype only if they share the block; a phased GT with no
    PS is treated as a single block, which is what chromosome-level
    phasing means. Multi-sample VCFs use the first sample column, the
    convention for a single-sample call set.
    """
    if len(fields) < 10 or not fields[8]:
        return "", frozenset()
    spec = dict(zip(fields[8].split(":"), fields[9].split(":")))
    gt = spec.get("GT", "")
    if "|" not in gt:
        return "", frozenset()
    haps = frozenset(i for i, a in enumerate(gt.split("|"))
                     if a == str(alt_index))
    if not haps:
        return "", frozenset()
    return spec.get("PS", "") or "*", haps


def parse_vcf(path: str | Path,
              require_pass: bool = True,
              keep_filters: tuple[str, ...] = ("PASS", "HighConf"),
              logger=None) -> Iterator[dict[str, Any]]:
    """Yield one record per ALT allele of a VCF.

    SEQC2 truth VCFs carry confidence in FILTER/INFO rather than QUAL, so
    `keep_filters` accepts both PASS and the tier labels used by that call
    set. Multi-allelic rows are split; symbolic ALTs (<DEL>, <DUP>) and
    breakend rows are skipped here because they are handled by the
    structural-variant / fusion path, not by codon substitution.
    """
    n_in = n_out = n_skip_symbolic = n_skip_filter = 0
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            n_in += 1
            chrom, pos, vid, ref, alts = f[0], int(f[1]), f[2], f[3].upper(), f[4]
            flt = f[6] if len(f) > 6 else "."
            info = f[7] if len(f) > 7 else ""
            if require_pass and flt not in (".", ""):
                labels = set(flt.split(";"))
                if not labels & set(keep_filters):
                    n_skip_filter += 1
                    continue
            for alt_i, alt in enumerate(alts.split(","), start=1):
                alt = alt.upper()
                if alt in (".", "*") or alt.startswith("<") or "[" in alt or "]" in alt:
                    n_skip_symbolic += 1
                    continue
                if not set(alt) <= set("ACGTN"):
                    n_skip_symbolic += 1
                    continue
                vclass = classify_small_variant(ref, alt)
                phase_set, haps = _phase(f, alt_i)
                n_out += 1
                yield {
                    "variant_id": f"{vclass}|{chrom}:{pos}{ref}>{alt}",
                    "variant_class": vclass,
                    "source": "sSNV_sIndel_VCF",
                    "genes": "",
                    "locus": f"{chrom}:{pos}",
                    "confidence": flt,
                    "payload": {
                        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                        "vcf_id": vid, "filter": flt, "info": info,
                        "phase_set": phase_set, "haplotypes": sorted(haps),
                        "af": _info_af(info),
                    },
                }
    if logger:
        logger.info(
            "VCF %s: %d rows -> %d alleles (%d filtered out, %d symbolic/skipped)",
            Path(path).name, n_in, n_out, n_skip_filter, n_skip_symbolic,
        )


# --------------------------------------------------------------------------
# 2. RNA editing sites (ANNOVAR *_multianno.txt)
# --------------------------------------------------------------------------

_AACHANGE_FIELDS = ("gene", "transcript", "exon", "cdna", "protein")


def parse_aachange(field: str) -> list[dict[str, str]]:
    """Split an ANNOVAR AAChange.refGene value into per-transcript records.

    Format: GENE:NM_x:exonN:c.A1525G:p.K509E , multiple entries comma-separated.
    Entries with fewer than 5 colon-separated parts (UTR/ncRNA rows) are
    returned with the fields that are present so nothing is silently lost.
    """
    out: list[dict[str, str]] = []
    if not field or field in (".", "UNKNOWN", "unknown"):
        return out
    for entry in field.split(","):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        rec = dict(zip(_AACHANGE_FIELDS, parts))
        rec["raw"] = entry
        out.append(rec)
    return out


def parse_annovar_res(path: str | Path, logger=None) -> Iterator[dict[str, Any]]:
    """Yield one record per RNA-editing site.

    The SEQC2 RES table repeats the VCF column names as a second header
    row (CHROM/POS/POS/REF/ALT); that row is detected and dropped.
    """
    n_in = n_out = 0
    with _open_text(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            n_in += 1
            chrom = (row.get("Chr") or "").strip()
            if not chrom or chrom.upper() == "CHROM":
                continue
            try:
                pos = int(row["Start"])
            except (KeyError, TypeError, ValueError):
                continue
            ref = (row.get("Ref") or "").upper()
            alt = (row.get("Alt") or "").upper()
            func = (row.get("Func.refGene") or ".").strip()
            exfunc = (row.get("ExonicFunc.refGene") or ".").strip()
            gene = (row.get("Gene.refGene") or "").strip()
            aac = parse_aachange(row.get("AAChange.refGene", ""))
            aac_ens = parse_aachange(row.get("AAChange.ensGene", ""))
            n_out += 1
            yield {
                "variant_id": f"RES|{chrom}:{pos}{ref}>{alt}",
                "variant_class": "RNA_EDITING",
                "source": "RES_ANNOVAR",
                "genes": gene,
                "locus": f"{chrom}:{pos}",
                "confidence": "",
                "payload": {
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "func_refgene": func,
                    "exonic_func_refgene": exfunc,
                    "gene_refgene": gene,
                    "gene_ensgene": (row.get("Gene.ensGene") or "").strip(),
                    "aachange_refgene": aac,
                    "aachange_ensgene": aac_ens,
                    "cosmic_coding": (row.get("cosmic95_coding") or ".").strip(),
                    "clnsig": (row.get("CLNSIG") or ".").strip(),
                    "interpro_domain": (row.get("Interpro_domain") or ".").strip(),
                },
            }
    if logger:
        logger.info("ANNOVAR RES %s: %d rows -> %d sites",
                    Path(path).name, n_in, n_out)


# --------------------------------------------------------------------------
# 3. Fusion genes
# --------------------------------------------------------------------------

_BP_RE = re.compile(r"^(chr[\w]+|[\w]+):(\d+)$")


def parse_fusions(path: str | Path, logger=None) -> Iterator[dict[str, Any]]:
    n_in = n_out = 0
    with _open_text(path) as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            n_in += 1
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            g1, g2 = row.get("Gene1", ""), row.get("Gene2", "")
            if not g1 or not g2:
                continue
            m1 = _BP_RE.match(row.get("breakpoint1", ""))
            m2 = _BP_RE.match(row.get("breakpoint2", ""))
            if not (m1 and m2):
                if logger:
                    logger.warning("fusion %s:%s has unparseable breakpoints "
                                   "(%r, %r) - skipped", g1, g2,
                                   row.get("breakpoint1"), row.get("breakpoint2"))
                continue
            n_out += 1
            try:
                conf = float(row.get("confidence", "") or "nan")
            except ValueError:
                conf = float("nan")
            yield {
                "variant_id": f"FUSION|{g1}--{g2}|{row['breakpoint1']}|{row['breakpoint2']}",
                "variant_class": "FUSION",
                "source": "Fusion_genes",
                "genes": f"{g1};{g2}",
                "locus": f"{row['breakpoint1']}::{row['breakpoint2']}",
                "confidence": conf,
                "payload": {
                    "gene1": g1, "gene2": g2,
                    "chrom1": m1.group(1), "pos1": int(m1.group(2)),
                    "chrom2": m2.group(1), "pos2": int(m2.group(2)),
                    "detection_type": row.get("Type", ""),
                    "validated": row.get("Validated", ""),
                    "tag": row.get("tag", ""),
                },
            }
    if logger:
        logger.info("fusions %s: %d rows -> %d events",
                    Path(path).name, n_in, n_out)


# --------------------------------------------------------------------------
# 4. Alternative splicing (SUPPA2-style event ids)
# --------------------------------------------------------------------------

AS_TYPES = {
    "SE": "skipping exon",
    "RI": "retained intron",
    "A3": "alternative 3' splice site",
    "A5": "alternative 5' splice site",
    "MX": "mutually exclusive exons",
    "AF": "alternative first exon",
    "AL": "alternative last exon",
}


def parse_suppa_event_id(event: str) -> dict[str, Any] | None:
    """Parse '<gene>;<TYPE>:<chr>:<coords...>:<strand>'.

    SUPPA2 encodes a variable number of coordinate groups depending on
    event type, so the groups are kept as a list and interpreted by the
    isoform builder rather than being force-fitted here.
    """
    if ";" not in event:
        return None
    gene, rest = event.split(";", 1)
    parts = rest.split(":")
    if len(parts) < 4:
        return None
    etype, chrom, strand = parts[0], parts[1], parts[-1]
    if strand not in ("+", "-"):
        return None
    groups: list[list[int]] = []
    for g in parts[2:-1]:
        try:
            groups.append([int(x) for x in g.split("-")])
        except ValueError:
            return None
    return {"gene_id": gene, "event_type": etype, "chrom": chrom,
            "strand": strand, "coord_groups": groups}


def parse_as_events(path: str | Path, logger=None) -> Iterator[dict[str, Any]]:
    n_in = n_out = n_bad = 0
    with _open_text(path) as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            n_in += 1
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            ev = row.get("AS", "")
            parsed = parse_suppa_event_id(ev)
            if parsed is None:
                n_bad += 1
                if logger:
                    logger.warning("unparseable AS event id: %r", ev)
                continue
            n_out += 1

            def _f(key: str) -> float:
                try:
                    return float(row.get(key, "") or "nan")
                except ValueError:
                    return float("nan")

            yield {
                "variant_id": f"AS_{parsed['event_type']}|{ev}",
                "variant_class": f"AS_{parsed['event_type']}",
                "source": "AS_LR",
                "genes": parsed["gene_id"],
                "locus": f"{parsed['chrom']}:{parsed['coord_groups'][0][0]}",
                "confidence": _f("confidence"),
                "payload": {
                    **parsed,
                    "event_id": ev,
                    "dPSI": _f("dPSI"),
                    "p_val": _f("p_val"),
                    "batch_support": row.get("Batch_support", ""),
                    "lib_support": row.get("Lib_support", ""),
                },
            }
    if logger:
        logger.info("AS events %s: %d rows -> %d events (%d unparseable)",
                    Path(path).name, n_in, n_out, n_bad)

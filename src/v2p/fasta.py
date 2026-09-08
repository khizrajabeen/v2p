"""FASTA emission with interchangeable header conventions.

The target header specification was not supplied with the input package,
so the writer is built around a registry of formatters. Switching
convention is a single `--header-style` flag and does not touch the
sequence-generation code, so the combined database can be re-emitted in
whatever layout the recipient specifies without re-running translation.

Every style carries the variant-type information that the request asked
for; they differ only in how it is encoded.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .build.smallvar import ProteinRecord
from .species import HUMAN, Species

# variant-type vocabulary emitted in headers (stable, documented values)
VARIANT_TYPE_VOCAB = {
    "SNV": "somatic single-nucleotide variant",
    "MNV": "somatic multi-nucleotide variant",
    "INDEL": "somatic insertion/deletion",
    "RNA_EDITING": "RNA editing site (A-to-I)",
    "FUSION": "gene fusion",
    "AS_SE": "alternative splicing: skipped exon",
    "AS_RI": "alternative splicing: retained intron",
    "AS_A3": "alternative splicing: alternative 3' splice site",
    "AS_A5": "alternative splicing: alternative 5' splice site",
    "AS_MX": "alternative splicing: mutually exclusive exons",
    "AS_AF": "alternative splicing: alternative first exon",
    "AS_AL": "alternative splicing: alternative last exon",
    "REFERENCE": "unmodified reference protein",
}

_ID_SAFE = re.compile(r"[^A-Za-z0-9_.:\-]+")


def sanitize_id(s: str) -> str:
    """Make an identifier safe for FASTA headers and search engines."""
    return _ID_SAFE.sub("_", s).strip("_") or "NA"


def wrap(seq: str, width: int = 60) -> str:
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# --------------------------------------------------------------------------
# header formatters
# --------------------------------------------------------------------------

def header_descriptive(r: ProteinRecord, prefix: str,
                       species: Species = HUMAN) -> str:
    """Pipe-delimited, self-documenting. Safe default."""
    fields = [
        f"{prefix}|{sanitize_id(r.seq_id)}",
        f"VT={r.variant_class}",
        f"CSQ={r.consequence or 'NA'}",
        f"GN={r.gene or 'NA'}",
        f"TX={r.transcript or 'NA'}",
        f"PC={r.protein_change or 'NA'}",
        f"LOC={r.locus or 'NA'}",
        f"SRC={r.source or 'NA'}",
        f"CONF={r.confidence or 'NA'}",
        f"LEN={len(r.sequence)}",
    ]
    if r.variant_pos_aa:
        fields.append(f"POS={r.variant_pos_aa}")
    if r.novel_span:
        fields.append(f"NOVEL={r.novel_span[0]}-{r.novel_span[1]}")
    if r.notes:
        fields.append("NOTE=" + ",".join(sanitize_id(n) for n in r.notes))
    return ">" + " ".join(fields)


def header_peff(r: ProteinRecord, prefix: str,
                species: Species = HUMAN) -> str:
    """PSI Extended FASTA Format entry line.

    Simple substitutions are encoded with \\VariantSimple=(pos|newAA);
    everything else (frameshift, fusion junction, isoform) has no
    residue-level PEFF representation, so the type is carried in the
    custom \\VariantType and \\Comment keys declared in the file header.
    """
    uid = sanitize_id(r.seq_id)
    parts = [f">{prefix}:{uid}", f"\\DbUniqueId={uid}"]
    pname = f"{r.gene or 'NA'} {VARIANT_TYPE_VOCAB.get(r.variant_class, r.variant_class)}"
    parts.append(f"\\PName={pname}")
    if r.gene:
        parts.append(f"\\GName={r.gene}")
    parts.append(f"\\TaxName={species.scientific_name}")
    parts.append(f"\\NcbiTaxId={species.taxon_id}")
    parts.append(f"\\Length={len(r.sequence)}")
    if r.consequence in ("missense", "stop_gained") and r.variant_pos_aa:
        new_aa = r.sequence[r.variant_pos_aa - 1] if r.variant_pos_aa <= len(r.sequence) else "X"
        parts.append(f"\\VariantSimple=({r.variant_pos_aa}|{new_aa})")
    parts.append(f"\\VariantType={r.variant_class}")
    parts.append(f"\\Consequence={r.consequence or 'NA'}")
    if r.transcript:
        parts.append(f"\\TranscriptId={r.transcript}")
    if r.protein_change:
        parts.append(f"\\ProteinChange={r.protein_change}")
    if r.locus:
        parts.append(f"\\GenomicLocus={r.locus}")
    if r.novel_span:
        parts.append(f"\\NovelSpan={r.novel_span[0]}-{r.novel_span[1]}")
    if r.source:
        parts.append(f"\\EvidenceSource={r.source}")
    if r.confidence:
        parts.append(f"\\Confidence={r.confidence}")
    if r.notes:
        parts.append("\\Comment=" + ";".join(r.notes))
    return " ".join(parts)


def header_pvac(r: ProteinRecord, prefix: str,
                species: Species = HUMAN) -> str:
    """Compact neoantigen-pipeline style: >MT.GENE.TX.consequence.change"""
    tag = "WT" if r.variant_class == "REFERENCE" else "MT"
    bits = [tag, r.gene or "NA", r.transcript or "NA",
            r.variant_class, r.consequence or "NA"]
    if r.protein_change:
        bits.append(r.protein_change.replace("p.", ""))
    return ">" + ".".join(sanitize_id(b) for b in bits)


def header_uniprot(r: ProteinRecord, prefix: str,
                   species: Species = HUMAN) -> str:
    """UniProt/SwissProt grammar, matching the supplied reference FASTA.

        >db|ACCESSION|ENTRY_NAME Description OS=Homo sapiens OX=9606 GN=SYM

    Reference entries keep `sp|` and their original accession so they stay
    byte-identical to the template. Variant entries use `vr|` with a
    generated accession, so a search engine can separate variant from
    canonical hits on the database prefix alone. Variant metadata is
    appended as further `KEY=value` pairs, which is exactly how UniProt
    already extends the line (OS=, OX=, GN=, PE=, SV=), so parsers that
    split on `\\b[A-Z]{2}=` keep working.
    """
    if r.variant_class == "REFERENCE":
        # A wild-type sequence is not a variant: emit it in plain SwissProt
        # form so it is indistinguishable from the supplied template, and
        # do not attach the triggering variant's locus to it.
        if r.extra.get("uniprot_header"):
            return r.extra["uniprot_header"]
        acc = r.extra.get("accession") or sanitize_id(r.seq_id)
        desc = r.extra.get("description") or f"{r.gene or 'NA'} reference protein"
        head = (f">sp|{acc}|"
                f"{species.entry_name(sanitize_id(r.gene or 'NA').upper())} "
                f"{desc} {species.os_ox()}")
        if r.gene:
            head += f" GN={r.gene}"
        head += " VT=REFERENCE CSQ=reference"
        if r.transcript:
            head += f" TX={r.transcript}"
        return head

    # A UniProt accession is not unique across our records (one gene can
    # carry several variants), so the entry id is a generated serial and the
    # UniProt accession is carried separately in UP=.
    acc = r.extra.get("_serial") or sanitize_id(r.seq_id)
    gene = r.gene or "NA"
    entry = (f"{species.entry_name(sanitize_id(gene).upper())}"
             f"_{r.variant_class}")
    desc = r.extra.get("description") or f"{gene} variant protein"
    change = r.protein_change or r.consequence or r.variant_class
    head = (f">vr|{prefix}_{acc}|{entry} {desc} ({change}; "
            f"{VARIANT_TYPE_VOCAB.get(r.variant_class, r.variant_class)})"
            f" {species.os_ox()}")
    if r.gene:
        head += f" GN={r.gene}"
    head += f" VT={r.variant_class} CSQ={r.consequence or 'NA'}"
    if r.extra.get("accession"):
        head += f" UP={r.extra['accession']}"
    if r.transcript:
        head += f" TX={r.transcript}"
    if r.protein_change:
        head += f" PC={r.protein_change}"
    if r.locus:
        head += f" LOC={r.locus}"
    if r.variant_pos_aa:
        head += f" POS={r.variant_pos_aa}"
    if r.novel_span:
        head += f" NOVEL={r.novel_span[0]}-{r.novel_span[1]}"
    if r.source:
        head += f" SRC={r.source}"
    if r.confidence:
        head += f" CONF={r.confidence}"
    if r.extra.get("uniprot_check"):
        head += f" QC={r.extra['uniprot_check']}"
    if r.extra.get("nmd") and r.extra["nmd"] != "not_applicable":
        head += f" NMD={r.extra['nmd']}"
    if r.extra.get("novel_peptides") is not None:
        head += f" NOVELPEP={r.extra['novel_peptides']}"
    if r.notes:
        head += " NOTE=" + ",".join(sanitize_id(n) for n in r.notes)
    return head


HEADER_STYLES: dict[str, Callable[..., str]] = {
    "descriptive": header_descriptive,
    "peff": header_peff,
    "pvac": header_pvac,
    "uniprot": header_uniprot,
}

# styles whose specification requires unwrapped sequence lines
SINGLE_LINE_STYLES = {"uniprot"}


# --------------------------------------------------------------------------
# writer
# --------------------------------------------------------------------------

def peff_file_header(db_name: str, prefix: str, n_entries: int,
                     description: str) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    lines = [
        "# PEFF 1.0",
        "# //",
        f"# DbName={db_name}",
        f"# DbDescription={description}",
        f"# Prefix={prefix}",
        "# DbSource=HCC1395 high-confidence variant package",
        f"# DbVersion={ts}",
        "# SequenceType=AA",
        f"# NumberOfEntries={n_entries}",
        "# GeneralComment=Custom keys VariantType, Consequence, TranscriptId, "
        "ProteinChange, GenomicLocus, NovelSpan, EvidenceSource, Confidence "
        "are non-standard extensions; PEFF-compliant readers ignore unknown keys.",
        "# //",
    ]
    return "\n".join(lines) + "\n"


def write_fasta(records: Sequence[ProteinRecord], path: str | Path,
                style: str = "descriptive", prefix: str = "HCC1395",
                width: int = 60, db_name: str = "HCC1395_variant_proteome",
                description: str = "HCC1395 variant protein database",
                dedup: bool = True, logger=None,
                species: Species = HUMAN) -> dict[str, int]:
    """Write records to FASTA. Returns per-variant-class counts.

    `species` defaults to human, so a caller that does not pass one emits
    exactly the headers it did before this parameter existed.
    """
    if style not in HEADER_STYLES:
        raise ValueError(f"unknown header style {style!r}; "
                         f"choose from {sorted(HEADER_STYLES)}")
    fmt = HEADER_STYLES[style]
    if style in SINGLE_LINE_STYLES:
        width = 0          # the reference template is unwrapped; match it
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    kept: list[ProteinRecord] = []
    seen_ids: set[str] = set()
    seen_seq: dict[str, ProteinRecord] = {}
    counts: dict[str, int] = {}
    merges: list[dict[str, str]] = []
    n_dup_seq = 0
    for r in records:
        if not r.sequence:
            continue
        base = sanitize_id(r.seq_id)
        uid = base
        i = 1
        while uid in seen_ids:
            i += 1
            uid = f"{base}_dup{i}"
        seen_ids.add(uid)
        r.seq_id = uid
        r.extra.setdefault(
            "_serial", f"{r.variant_class}_{len(seen_ids):06d}")
        if dedup:
            key = f"{r.variant_class}\t{r.sequence}"
            if key in seen_seq:
                n_dup_seq += 1
                # Dropping the record silently loses which variant it came
                # from. Two different DNA changes encoding the same residue
                # collapse here, and only the first one's locus survives
                # into the header - which made a benchmark score the second
                # as a miss when its protein had in fact been built.
                # Record the merge instead of discarding it.
                _k = seen_seq[key]
                merges.append({
                    "kept_seq_id": _k.seq_id,
                    # The join key stage 6 can actually use: the header
                    # carries a generated serial, not this seq_id, so
                    # locus+transcript is what the two sides share.
                    "kept_locus": _k.locus or "",
                    "kept_transcript": _k.transcript or "",
                    "merged_seq_id": uid,
                    "merged_locus": r.locus or "",
                    "merged_variant_class": r.variant_class,
                    "merged_protein_change": r.protein_change or "",
                    "merged_transcript": r.transcript or "",
                })
                continue
            seen_seq[key] = r
        kept.append(r)
        counts[r.variant_class] = counts.get(r.variant_class, 0) + 1

    with open(path, "w", encoding="utf-8") as fh:
        if style == "peff":
            fh.write(peff_file_header(db_name, prefix, len(kept), description))
        for r in kept:
            fh.write(fmt(r, prefix, species) + "\n")
            fh.write(wrap(r.sequence, width) + "\n")

    # Sidecar, not header bytes: the merge record must not change the
    # FASTA, which is compared byte-for-byte across releases.
    if merges:
        mpath = Path(str(path) + ".merged_loci.tsv")
        with open(mpath, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(merges[0]),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(merges)

    if logger:
        logger.info("wrote %d sequences to %s (style=%s, %d identical "
                    "sequences collapsed)", len(kept), path, style, n_dup_seq)
        for k in sorted(counts):
            logger.info("  fasta class %s = %d", k, counts[k])
    counts["_total"] = len(kept)
    counts["_dedup_collapsed"] = n_dup_seq
    return counts


def write_record_table(records: Iterable[ProteinRecord],
                       path: str | Path) -> int:
    """Companion TSV: one row per FASTA entry, for joins and QC."""
    cols = ["seq_id", "variant_class", "consequence", "gene", "transcript",
            "protein_change", "variant_pos_aa", "protein_length",
            "ref_protein_len", "locus", "source", "confidence",
            "novel_span", "notes"]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            n += 1
            fh.write("\t".join([
                r.seq_id, r.variant_class, r.consequence, r.gene, r.transcript,
                r.protein_change,
                str(r.variant_pos_aa) if r.variant_pos_aa else "",
                str(len(r.sequence)),
                str(r.ref_protein_len) if r.ref_protein_len is not None else "",
                r.locus, r.source, str(r.confidence),
                f"{r.novel_span[0]}-{r.novel_span[1]}" if r.novel_span else "",
                ";".join(r.notes),
            ]) + "\n")
    return n

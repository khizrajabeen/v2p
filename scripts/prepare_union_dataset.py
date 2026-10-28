#!/usr/bin/env python3
"""Convert a Variant_Union_Dataset tree into a v2p input folder.

    SNV/*.csv.gz + InDel/*.csv   ->  somatic_snv_indel.vcf
    RNA_editing/*.tsv            ->  rna_editing.txt
    Gene_Fusion_{SR,LR}/*.tsv    ->  fusions.csv
    AS/*.tsv                     ->  splicing.csv

None of the union tables is in a shape v2p detects (docs/USAGE.md step
2), so each is rewritten rather than reinterpreted at run time. Two
things are corrected on the way through:

  * Contig naming. The tables use UCSC scaffold names
    (chr14_GL000009v2_random); GENCODE's primary assembly calls the same
    scaffold GL000009.2. The accession is recoverable, so it is mapped.
    What is still absent from the reference index is dropped and counted.

  * Support. A union over every tool and file keeps events seen once.
    --min-tools / --min-batches select a consensus subset; the defaults
    keep everything.

SV/SV_exact_union.tsv is not converted: v2p consumes four evidence types
and a breakend pair is none of them. Breakends joining two genes are
already in the fusion tables.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

csv.field_size_limit(1 << 30)

# UCSC random/unplaced scaffold -> GENCODE (Ensembl) accession.
#   chr14_GL000009v2_random -> GL000009.2
#   chrUn_KI270742v1        -> KI270742.1
_UCSC_SCAFFOLD = re.compile(
    r"^chr(?:[\dXYM]+|Un)_([A-Z]+\d+)v(\d+)(?:_(?:random|alt|fix))?$")


def contig_mapper(fai: Path):
    """Return (map_fn, known_names) for the contigs of one reference."""
    known = {line.split("\t", 1)[0] for line in fai.read_text().splitlines()
             if line.strip()}

    def to_ref(chrom: str) -> str | None:
        c = chrom.strip()
        if c in known:
            return c
        m = _UCSC_SCAFFOLD.match(c)
        if m:
            cand = f"{m.group(1)}.{m.group(2)}"
            if cand in known:
                return cand
        for cand in (c[3:] if c.startswith("chr") else "chr" + c,):
            if cand in known:
                return cand
        if c in ("chrM", "chrMT", "MT", "M"):
            for cand in ("chrM", "MT", "M", "chrMT"):
                if cand in known:
                    return cand
        return None

    return to_ref, known


def _open_text(p: Path):
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace",
                         newline="")
    return open(p, "rt", encoding="utf-8", errors="replace", newline="")


def _int(row: dict, key: str, default: int = 0) -> int:
    try:
        return int(float(row.get(key) or default))
    except (TypeError, ValueError):
        return default


def _pct(part: int, whole: int) -> str:
    """A support fraction in [0,1], or '' when the denominator is unknown.

    `confidence` reaches the FASTA header verbatim, so '' is the honest
    value for "not stated"; 0.0 would read as "no support".
    """
    return f"{part / whole:.4f}" if whole > 0 else ""


# --------------------------------------------------------------------------
# small variants -> one VCF
# --------------------------------------------------------------------------

_CONTIG_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def convert_small_variants(snv: Path | None, indel: Path | None,
                           out: Path, to_ref, min_tools: int,
                           min_batches: int, genome=None) -> dict:
    stats = Counter()
    rows: list[tuple[int, int, str, str, str, str]] = []
    order = {c: i for i, c in enumerate(_CONTIG_ORDER)}

    for path, kind in ((snv, "SNV"), (indel, "INDEL")):
        if path is None:
            continue
        with _open_text(path) as fh:
            for row in csv.DictReader(fh):
                stats[f"{kind}_in"] += 1
                chrom = to_ref(row.get("CHROM", ""))
                if chrom is None:
                    stats[f"{kind}_dropped_contig"] += 1
                    continue
                tools = _int(row, "detected_tool_count")
                batches = _int(row, "detected_batch_count")
                if tools < min_tools or batches < min_batches:
                    stats[f"{kind}_dropped_support"] += 1
                    continue
                try:
                    pos = int(row["POS"])
                except (KeyError, ValueError):
                    stats[f"{kind}_dropped_badpos"] += 1
                    continue
                ref, alt = (row.get("REF") or "").upper(), (row.get("ALT") or "").upper()
                if not ref or not alt:
                    stats[f"{kind}_dropped_noallele"] += 1
                    continue
                if not (set(ref) <= set("ACGTN") and set(alt) <= set("ACGTN")):
                    stats[f"{kind}_dropped_nonacgt"] += 1
                    continue
                info = [
                    f"TOOLS={tools}",
                    f"FILES={_int(row, 'detected_file_count')}",
                    f"BATCHES={batches}",
                    f"TOOLSET={_toolset(row)}",
                ]
                fsr = (row.get("file_support_rate_percent") or "").strip()
                if fsr:
                    info.append(f"FSR={fsr}")
                bsr = (row.get("batch_support_rate_percent") or "").strip()
                if bsr:
                    info.append(f"BSR={bsr}")
                vtype = (row.get("variant_type") or "").strip()
                if vtype:
                    info.append(f"VTYPE={vtype}")
                rows.append((order.get(chrom, 99 + len(order)), pos,
                             chrom, ref, alt, ";".join(info)))
                stats[f"{kind}_out"] += 1

    rows.sort(key=lambda r: (r[0], r[1], r[3], r[4]))

    n_mismatch = 0
    if genome is not None:
        for _, pos, chrom, ref, _alt, _info in rows:
            got = genome.fetch(chrom, pos, pos + len(ref) - 1)
            if got and got.upper() != ref:
                n_mismatch += 1
        stats["ref_mismatch"] = n_mismatch

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##reference=GRCh38\n")
        fh.write("##source=Variant_Union_Dataset (SNV + InDel union tables)\n")
        if genome is not None:
            for c in _CONTIG_ORDER:
                ln = genome.length(c)
                if ln:
                    fh.write(f"##contig=<ID={c},length={ln}>\n")
        for key, num, typ, desc in (
            ("TOOLS", "1", "Integer", "callers detecting this variant"),
            ("FILES", "1", "Integer", "input files detecting this variant"),
            ("BATCHES", "1", "Integer", "sequencing batches detecting this variant"),
            ("TOOLSET", "1", "String", "names of the detecting callers"),
            ("FSR", "1", "Float", "file support rate, percent, from the union table"),
            ("BSR", "1", "Float", "batch support rate, percent, from the union table"),
            ("VTYPE", "1", "String", "variant_type as given by the union table"),
        ):
            fh.write(f'##INFO=<ID={key},Number={num},Type={typ},'
                     f'Description="{desc}">\n')
        fh.write("##FILTER=<ID=PASS,Description=\"in the union call set\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, pos, chrom, ref, alt, info in rows:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
    return stats


def _toolset(row: dict) -> str:
    """Which of the three callers supports this row."""
    names = []
    for tool in ("Mutect2", "Strelka2", "DeepSomatic"):
        if _int(row, f"{tool}_file_count") > 0:
            names.append(tool)
    return "|".join(names) or "."


# --------------------------------------------------------------------------
# RNA editing -> ANNOVAR-style table
# --------------------------------------------------------------------------

_RES_COLUMNS = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
                "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene",
                "AAChange.refGene", "edit_type", "file_support_count",
                "batch_support_count", "tool_support_pattern"]


def convert_rna_editing(src: Path, out: Path, to_ref, min_files: int,
                        min_batches: int) -> dict:
    """Rewrite the A-to-I union table in the shape v2p detects.

    Annotation columns are emitted as '.' rather than invented: v2p needs
    only Chr/Start/Ref/Alt and derives the consequence itself, while
    `Func.refGene` has to exist for the file to be recognised. Stage 3
    counts recoding sites from `AAChange.refGene`, so a guess there would
    be reported as a measurement.
    """
    stats = Counter()
    out.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(src) as fh, open(out, "w", encoding="utf-8",
                                     newline="\n") as fo:
        fo.write("\t".join(_RES_COLUMNS) + "\n")
        for row in csv.DictReader(fh, delimiter="\t"):
            stats["in"] += 1
            edit = (row.get("edit_type") or "").strip().upper()
            if edit not in ("A>G", "T>C"):
                stats["dropped_not_AtoI"] += 1
                continue
            chrom = to_ref(row.get("CHROM", ""))
            if chrom is None:
                stats["dropped_contig"] += 1
                continue
            files = _int(row, "file_support_count")
            batches = _int(row, "batch_support_count")
            if files < min_files or batches < min_batches:
                stats["dropped_support"] += 1
                continue
            try:
                pos = int(row["POS"])
            except (KeyError, ValueError):
                stats["dropped_badpos"] += 1
                continue
            ref = (row.get("REF") or "").upper()
            alt = (row.get("ALT") or "").upper()
            if len(ref) != 1 or len(alt) != 1 or ref == alt:
                stats["dropped_not_substitution"] += 1
                continue
            fo.write("\t".join([
                chrom, str(pos), str(pos), ref, alt,
                ".", ".", ".", ".", ".",
                edit, str(files), str(batches),
                (row.get("global_tool_support_pattern") or ".").strip() or ".",
            ]) + "\n")
            stats["out"] += 1
    return stats


# --------------------------------------------------------------------------
# fusions -> CSV with chr:pos breakpoints
# --------------------------------------------------------------------------

def convert_fusions(sr: Path | None, lr: Path | None, out: Path, to_ref,
                    min_tools: int, min_batches: int) -> dict:
    """Merge the short-read and long-read fusion unions into one table.

    Partner order is kept as the callers reported it; they all name the
    5' partner first. `Type` records which read length saw the event, and
    an event seen by both becomes one row tagged `SR;LR`.
    """
    stats = Counter()
    events: dict[tuple[str, str, str, int, str, int], dict] = {}
    # Taken over every row read, not the rows kept, so `confidence` does
    # not depend on the threshold this was run with.
    max_batches = 0

    for path, kind in ((sr, "SR"), (lr, "LR")):
        if path is None:
            continue
        with _open_text(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                stats[f"{kind}_in"] += 1
                g1 = (row.get("Gene1") or "").strip()
                g2 = (row.get("Gene2") or "").strip()
                if not g1 or not g2 or g1 == "." or g2 == ".":
                    stats[f"{kind}_dropped_nogene"] += 1
                    continue
                c1 = to_ref(row.get("Chrom1", ""))
                c2 = to_ref(row.get("Chrom2", ""))
                if c1 is None or c2 is None:
                    stats[f"{kind}_dropped_contig"] += 1
                    continue
                try:
                    p1 = int(row["Breakpoint1"])
                    p2 = int(row["Breakpoint2"])
                except (KeyError, ValueError):
                    stats[f"{kind}_dropped_badbp"] += 1
                    continue
                tools = _int(row, "detected_tool_count")
                batches = _int(row, "detected_batch_count")
                max_batches = max(max_batches, batches)
                if tools < min_tools or batches < min_batches:
                    stats[f"{kind}_dropped_support"] += 1
                    continue
                key = (g1, g2, c1, p1, c2, p2)
                prev = events.get(key)
                if prev is None:
                    events[key] = {
                        "Gene1": g1, "Gene2": g2,
                        "breakpoint1": f"{c1}:{p1}",
                        "breakpoint2": f"{c2}:{p2}",
                        "Type": kind,
                        "tools": tools, "batches": batches,
                        "files": _int(row, "detected_file_count"),
                        "toolnames": (row.get("detected_tools") or "").strip(),
                        "strand1": (row.get("Strand1_set") or "").strip(),
                        "strand2": (row.get("Strand2_set") or "").strip(),
                    }
                    stats[f"{kind}_new"] += 1
                else:
                    if kind not in prev["Type"].split(";"):
                        prev["Type"] = ";".join(sorted(
                            set(prev["Type"].split(";")) | {kind}))
                    prev["tools"] = max(prev["tools"], tools)
                    prev["batches"] = max(prev["batches"], batches)
                    prev["files"] = max(prev["files"], _int(row, "detected_file_count"))
                    stats[f"{kind}_merged"] += 1

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, lineterminator="\n")
        w.writerow(["tag", "Gene1", "Gene2", "Type", "Validated",
                    "confidence", "breakpoint1", "breakpoint2",
                    "tools", "batches", "files", "detected_tools",
                    "strand1", "strand2"])
        for key in sorted(events):
            e = events[key]
            w.writerow([
                f"{e['Gene1']}--{e['Gene2']}", e["Gene1"], e["Gene2"],
                e["Type"], "", _pct(e["batches"], max_batches),
                e["breakpoint1"], e["breakpoint2"],
                e["tools"], e["batches"], e["files"], e["toolnames"],
                e["strand1"], e["strand2"],
            ])
            stats["out"] += 1
    stats["max_batches"] = max_batches
    return stats


# --------------------------------------------------------------------------
# splicing -> CSV with an AS column
# --------------------------------------------------------------------------

def convert_splicing(src: Path, out: Path, to_ref, min_files: int,
                     min_batches: int) -> dict:
    """Rewrite the SUPPA2 event union as a CSV v2p can read.

    `dPSI` and `p_val` are left empty: this is a union of events observed
    in the tumour, not a differential test, so there is no delta-PSI.
    """
    stats = Counter()
    max_batches = 0
    kept: list[tuple[str, str, int, int, str]] = []
    with _open_text(src) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            stats["in"] += 1
            ev = (row.get("event_id") or "").strip()
            if ";" not in ev or ":" not in ev:
                stats["dropped_badevent"] += 1
                continue
            # '<gene>;<TYPE>:<chr>:<coords>:<strand>'
            body = ev.split(";", 1)[1].split(":")
            if len(body) < 4:
                stats["dropped_badevent"] += 1
                continue
            chrom = to_ref(body[1])
            if chrom is None:
                stats["dropped_contig"] += 1
                continue
            files = _int(row, "file_support_count")
            batches = _int(row, "batch_support_count")
            max_batches = max(max_batches, batches)
            if files < min_files or batches < min_batches:
                stats["dropped_support"] += 1
                continue
            kept.append((ev, (row.get("gene_id") or "").strip(), files,
                         batches, (row.get("tool_support_pattern") or "").strip()))

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, lineterminator="\n")
        w.writerow(["AS", "dPSI", "p_val", "Batch_support", "Lib_support",
                    "gene_id", "confidence", "tool_support_pattern"])
        for ev, gid, files, batches, pattern in kept:
            w.writerow([ev, "", "", batches, files, gid,
                        _pct(batches, max_batches), pattern])
            stats["out"] += 1
    stats["max_batches"] = max_batches
    return stats


# --------------------------------------------------------------------------

class _Ref:
    """Minimal pyfaidx front end: contig lengths and REF-allele checks."""

    def __init__(self, fasta: Path):
        from pyfaidx import Fasta
        self.fa = Fasta(str(fasta), sequence_always_upper=True, as_raw=True)
        self.names = set(self.fa.keys())

    def length(self, chrom: str) -> int | None:
        return len(self.fa[chrom]) if chrom in self.names else None

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.names:
            return ""
        try:
            return str(self.fa[chrom][start - 1:end])
        except Exception:
            return ""


def find(root: Path, *patterns: str) -> Path | None:
    for pat in patterns:
        hits = sorted(root.glob(pat))
        if hits:
            return hits[0]
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", type=Path,
                    help="the extracted Variant_Union_Dataset directory")
    ap.add_argument("--outdir", type=Path, required=True,
                    help="v2p input folder to write")
    ap.add_argument("--fai", type=Path, required=True,
                    help="reference .fai, for contig names and order")
    ap.add_argument("--genome", type=Path,
                    help="reference FASTA; enables ##contig headers and a "
                         "REF-allele check on the small variants")
    ap.add_argument("--min-tools", type=int, default=1,
                    help="keep small variants / fusions seen by at least "
                         "this many callers (default 1, i.e. keep all)")
    ap.add_argument("--min-batches", type=int, default=1,
                    help="keep events seen in at least this many sequencing "
                         "batches (default 1, i.e. keep all)")
    ap.add_argument("--min-files", type=int, default=1,
                    help="keep RNA-editing sites and splicing events seen in "
                         "at least this many files (default 1)")
    ap.add_argument("--skip", action="append", default=[],
                    choices=["small", "editing", "fusion", "splicing"],
                    help="skip one converter; repeatable")
    args = ap.parse_args(argv)

    ds: Path = args.dataset
    if not ds.is_dir():
        print(f"not a directory: {ds}", file=sys.stderr)
        return 2
    to_ref, known = contig_mapper(args.fai)
    genome = _Ref(args.genome) if args.genome else None
    print(f"reference declares {len(known)} contigs")

    report: dict[str, Counter] = {}

    if "small" not in args.skip:
        snv = find(ds, "SNV/*.csv.gz", "SNV/*.csv")
        indel = find(ds, "InDel/*.csv", "INDEL/*.csv", "Indel/*.csv")
        if snv or indel:
            print("converting small variants ...", flush=True)
            report["small_variants"] = convert_small_variants(
                snv, indel, args.outdir / "somatic_snv_indel.vcf", to_ref,
                args.min_tools, args.min_batches, genome)

    if "editing" not in args.skip:
        res = find(ds, "RNA_editing/*.tsv", "RNA_Editing/*.tsv")
        if res:
            print("converting RNA editing ...", flush=True)
            report["rna_editing"] = convert_rna_editing(
                res, args.outdir / "rna_editing.txt", to_ref,
                args.min_files, args.min_batches)

    if "fusion" not in args.skip:
        sr = find(ds, "Gene_Fusion_SR/*.tsv")
        lr = find(ds, "Gene_Fusion_LR/*.tsv")
        if sr or lr:
            print("converting fusions ...", flush=True)
            report["fusions"] = convert_fusions(
                sr, lr, args.outdir / "fusions.csv", to_ref,
                args.min_tools, args.min_batches)

    if "splicing" not in args.skip:
        as_ = find(ds, "AS/*.tsv")
        if as_:
            print("converting splicing ...", flush=True)
            report["splicing"] = convert_splicing(
                as_, args.outdir / "splicing.csv", to_ref,
                args.min_files, args.min_batches)

    print()
    for name, st in report.items():
        print(f"== {name}")
        for key in sorted(st):
            print(f"   {key:28s} {st[key]:>12,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

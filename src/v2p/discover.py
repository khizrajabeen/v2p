"""Work out what is in an input folder, without relying on filenames.

The same data arrives named differently from every collaborator, so every
file is identified by its content. Each decision carries the evidence
that produced it, and ties are reported rather than guessed at: a wrong
guess here produces a plausible database from the wrong data, which is
the worst failure mode this pipeline has.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

# evidence types the pipeline can consume
SMALL_VARIANTS = "small_variants"
RNA_EDITING = "rna_editing"
FUSION = "fusion"
SPLICING = "splicing"
REGIONS = "regions"
GENOME = "genome"
ANNOTATION = "annotation"
PROTEOME = "proteome"
TRANSLATIONS = "translations"
INDEX = "index"
UNKNOWN = "unknown"

ROLE_LABEL = {
    SMALL_VARIANTS: "somatic SNV / InDel calls",
    RNA_EDITING: "RNA editing sites",
    FUSION: "gene fusion calls",
    SPLICING: "alternative splicing events",
    REGIONS: "high-confidence regions",
    GENOME: "reference genome",
    ANNOTATION: "gene annotation",
    PROTEOME: "reference proteome",
    TRANSLATIONS: "annotation-source protein translations",
    INDEX: "index file (no data of its own)",
    UNKNOWN: "unrecognised",
}

_SUPPA = re.compile(r"^[^;]+;(SE|RI|MX|A3|A5|AF|AL):[^:]+:.+:[+-]$")
_LOCUS = re.compile(r"^(chr)?[\w.]+:\d+$")


@dataclass
class Detection:
    path: Path
    role: str
    confidence: float
    evidence: list[str] = field(default_factory=list)
    n_records: int | None = None
    detail: dict = field(default_factory=dict)

    def __str__(self) -> str:
        n = f", {self.n_records:,} records" if self.n_records else ""
        return (f"{self.path.name}: {ROLE_LABEL[self.role]} "
                f"({self.confidence:.0%}{n})")


def _open(p: Path):
    with open(p, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "rt", encoding="utf-8", errors="replace")


def _head(p: Path, n_bytes: int = 1 << 18) -> str:
    try:
        with _open(p) as fh:
            return fh.read(n_bytes)
    except Exception:
        return ""


def _sniff_delimiter(line: str) -> str:
    return "\t" if line.count("\t") > line.count(",") else ","


def _count_records(p: Path, is_vcf: bool = False, cap: int = 5_000_000) -> int:
    n = 0
    try:
        with _open(p) as fh:
            for line in fh:
                if is_vcf and line.startswith("#"):
                    continue
                n += 1
                if n >= cap:
                    break
    except Exception:
        return 0
    return n if is_vcf else max(0, n - 1)


def detect_file(path: Path) -> Detection:
    """Identify one file by content."""
    p = Path(path)
    ev: list[str] = []

    with open(p, "rb") as fh:
        magic = fh.read(4)
    if magic[:4] in (b"CSI\x01", b"TBI\x01"):
        return Detection(p, INDEX, 1.0, ["CSI/TBI magic bytes"])
    if magic[:2] == b"PK":
        return Detection(p, UNKNOWN, 1.0, ["zip archive - extract it first"])

    # Companion index files sit next to the data they index and often
    # share its shape: a .fai is tab-delimited with numeric columns 2 and
    # 3, which is exactly a BED. Identify them by their relationship to a
    # neighbouring file rather than by their contents.
    if p.suffix in (".fai", ".gzi", ".bai", ".crai", ".tbi", ".csi", ".idx"):
        base = p.with_suffix("")
        rel = (f"index for {base.name}" if base.exists()
               else f"{p.suffix} index file")
        return Detection(p, INDEX, 1.0, [rel])

    head = _head(p)
    # CSI and TBI indexes are BGZF-compressed, so their own magic bytes sit
    # behind the gzip header and only appear after decompression.
    if head[:4] in ("CSI\x01", "TBI\x01"):
        return Detection(p, INDEX, 1.0,
                         ["CSI/TBI magic bytes inside a BGZF container"])
    if not head:
        return Detection(p, UNKNOWN, 0.0, ["unreadable or binary"])

    # ---- VCF -----------------------------------------------------------
    if head.startswith("##fileformat=VCF"):
        ev.append("##fileformat=VCF header")
        det = {}
        m = re.search(r"##reference=(\S+)", head)
        if m:
            det["reference"] = m.group(1)
            ev.append(f"declares reference {m.group(1)}")
        contigs = re.findall(r"##contig=<ID=([^,>]+),length=(\d+)", head)
        if contigs:
            det["contigs"] = len(contigs)
            for c, L in contigs:
                if c in ("chr1", "1"):
                    det["build"] = ("GRCh38" if L == "248956422"
                                    else "GRCh37" if L == "249250621"
                                    else f"unknown (chr1={L})")
                    ev.append(f"chr1 length implies {det['build']}")
                    break
        return Detection(p, SMALL_VARIANTS, 1.0, ev,
                         _count_records(p, is_vcf=True), det)

    # ---- FASTA ---------------------------------------------------------
    if head.lstrip().startswith(">"):
        body = "".join(l for l in head.splitlines()[1:40]
                       if not l.startswith(">"))
        letters = set(body.upper()) - set("\n\r ")
        nucleotide = letters and letters <= set("ACGTUN")
        n_entries = head.count("\n>") + (1 if head.startswith(">") else 0)
        if nucleotide:
            ev.append("nucleotide alphabet")
            return Detection(p, GENOME, 0.95, ev, None,
                             {"first_contig": head[1:head.find("\n")][:40]})
        ev.append("protein alphabet")
        if re.search(r"^>\w+\|\w+\|\w+_\w+ .*\bOS=", head, re.M):
            ev.append("UniProt header grammar (db|accession|entry ... OS=)")
            return Detection(p, PROTEOME, 0.98, ev, _count_fasta(p))
        if re.search(r"^>ENSP\d+\.\d+\|ENST", head, re.M):
            ev.append("ENSP|ENST pipe-delimited header")
            return Detection(p, TRANSLATIONS, 0.98, ev, _count_fasta(p))
        return Detection(p, PROTEOME, 0.6,
                         ev + ["protein FASTA, header style unrecognised"],
                         _count_fasta(p))

    # ---- GTF / GFF -----------------------------------------------------
    if 'transcript_id "' in head or "transcript_id=" in head:
        ev.append("transcript_id attribute")
        feats = {f.split("\t")[2] for f in head.splitlines()
                 if not f.startswith("#") and f.count("\t") >= 8}
        return Detection(p, ANNOTATION, 0.97, ev, None,
                         {"features": sorted(feats)[:8]})

    # ---- delimited text ------------------------------------------------
    first = head.splitlines()[0] if head.splitlines() else ""
    delim = _sniff_delimiter(first)
    try:
        rdr = csv.reader(io.StringIO(head), delimiter=delim)
        rows = [r for _, r in zip(range(60), rdr)]
    except Exception:
        rows = []
    if not rows:
        return Detection(p, UNKNOWN, 0.0, ["no parseable structure"])

    header = [h.strip() for h in rows[0]]
    lower = [h.lower() for h in header]
    data = rows[1:]

    # ANNOVAR RNA-editing table
    if any(h.startswith("Func.") for h in header) and "Chr" in header:
        ev.append("ANNOVAR annotation columns (Func.*, Chr, Start, Ref, Alt)")
        alleles = set()
        try:
            i = {h: k for k, h in enumerate(header)}
            for r in data[:200]:
                if len(r) <= max(i["Ref"], i["Alt"]):
                    continue
                ref_a, alt_a = r[i["Ref"]].strip(), r[i["Alt"]].strip()
                # ANNOVAR output from some pipelines repeats the source VCF
                # column names as a second header row; skip anything that is
                # not an actual allele rather than counting it as one.
                if not ref_a or not alt_a:
                    continue
                if not set(ref_a.upper()) <= set("ACGTN-"):
                    continue
                if not set(alt_a.upper()) <= set("ACGTN-"):
                    continue
                alleles.add(f"{ref_a.upper()}>{alt_a.upper()}")
        except KeyError:
            pass
        if alleles and alleles <= {"A>G", "T>C"}:
            ev.append("all sampled alleles are A>G or T>C, as A-to-I editing predicts")
            return Detection(p, RNA_EDITING, 0.98, ev, _count_records(p))
        return Detection(p, RNA_EDITING, 0.75,
                         ev + ["allele spectrum not exclusively A>G/T>C"],
                         _count_records(p))

    # SUPPA2 splicing events - look for the event-id grammar in any column
    for col in range(min(4, len(header))):
        vals = [r[col] for r in data[:40] if len(r) > col]
        if vals and sum(bool(_SUPPA.match(v.strip())) for v in vals) >= max(2, len(vals) // 2):
            ev.append(f"column '{header[col]}' matches the SUPPA2 event grammar")
            return Detection(p, SPLICING, 0.97, ev, _count_records(p),
                             {"event_column": header[col]})

    # fusion table - two gene columns and two chr:pos columns
    gene_cols = [h for h in lower if re.fullmatch(r"gene ?[12]|gene_?[ab]|"
                                                 r"[53]'?_?(prime_?)?gene", h)]
    locus_cols = []
    for col, name in enumerate(header):
        vals = [r[col].strip() for r in data[:30] if len(r) > col]
        if vals and sum(bool(_LOCUS.match(v)) for v in vals) >= max(2, len(vals) // 2):
            locus_cols.append(name)
    if len(gene_cols) >= 2 and len(locus_cols) >= 2:
        ev.append(f"gene columns {gene_cols} and locus columns {locus_cols}")
        return Detection(p, FUSION, 0.95, ev, _count_records(p),
                         {"gene_columns": gene_cols, "locus_columns": locus_cols})

    # BED
    if delim == "\t" and len(header) >= 3:
        ok = 0
        for r in rows[:30]:
            if len(r) >= 3 and r[1].isdigit() and r[2].isdigit():
                ok += 1
        if ok >= max(3, len(rows[:30]) - 2):
            ev.append("three or more columns, second and third integer")
            return Detection(p, REGIONS, 0.9, ev, _count_records(p) + 1)

    return Detection(p, UNKNOWN, 0.0,
                     [f"delimited ({delim!r}) with columns {header[:6]}"])


def _count_fasta(p: Path) -> int:
    n = 0
    try:
        with _open(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
    except Exception:
        return 0
    return n


def detect_folder(folder: str | Path, recursive: bool = True) -> list[Detection]:
    """Identify every file in a folder. Deterministic order."""
    root = Path(folder)
    files = sorted(f for f in (root.rglob("*") if recursive else root.iterdir())
                   if f.is_file() and not f.name.startswith("."))
    return [detect_file(f) for f in files]


def build_run_plan(dets: list[Detection]) -> dict:
    """Turn detections into a run plan, flagging ambiguity.

    Where several files claim the same role the highest confidence wins
    and the rest are listed as alternatives, because silently picking one
    is how the wrong VCF ends up in a release.
    """
    by_role: dict[str, list[Detection]] = {}
    for d in dets:
        if d.role in (UNKNOWN, INDEX):
            continue
        by_role.setdefault(d.role, []).append(d)

    def _stem(path: Path) -> str:
        """Identity ignoring compression: X and X.gz are the same file."""
        n = path.name
        for suf in (".gz", ".bgz", ".bz2", ".zst"):
            if n.endswith(suf):
                n = n[: -len(suf)]
                break
        return n

    plan: dict = {"selected": {}, "alternatives": {}, "problems": []}
    for role, ds in by_role.items():
        # Prefer uncompressed, then confidence, then record count. Random
        # access into a plain FASTA is far cheaper than into a gzip, and
        # the two hold identical data.
        ds.sort(key=lambda d: (d.path.suffix in (".gz", ".bgz", ".bz2", ".zst"),
                               -d.confidence, -(d.n_records or 0), d.path.name))
        plan["selected"][role] = str(ds[0].path)
        if len(ds) > 1:
            plan["alternatives"][role] = [str(x.path) for x in ds[1:]]
            # Only genuinely different files are ambiguous. The same file
            # in two compression states is not a decision the user needs
            # to make, and asking makes real warnings easier to ignore.
            distinct = {_stem(d.path) for d in ds}
            if len(distinct) > 1 and ds[0].confidence - ds[1].confidence < 0.05:
                plan["problems"].append(
                    f"{ROLE_LABEL[role]}: {len(distinct)} candidates of "
                    f"similar confidence - name the one you want explicitly")

    evidence = [r for r in (SMALL_VARIANTS, RNA_EDITING, FUSION, SPLICING)
                if r in plan["selected"]]
    if not evidence:
        plan["problems"].append("no variant evidence found in this folder")
    plan["evidence_types"] = evidence

    for role, label in ((GENOME, "reference genome"),
                        (ANNOTATION, "gene annotation")):
        if role not in plan["selected"]:
            plan["problems"].append(f"no {label} found - supply it with --ref")

    # build consistency, where it can be determined
    builds = {d.detail.get("build") for d in dets if d.detail.get("build")}
    builds.discard(None)
    if len(builds) > 1:
        plan["problems"].append(f"inputs declare different genome builds: {builds}")
    elif builds:
        plan["genome_build"] = builds.pop()

    return plan

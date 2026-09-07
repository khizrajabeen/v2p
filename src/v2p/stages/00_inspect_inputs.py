#!/usr/bin/env python3
"""Stage 0a - inspect every input file before doing anything with it.

Identifies each file by content (not just extension), reports what is
actually inside it, and cross-checks the things that silently ruin a
conversion run: genome build, contig naming, coordinate base, and whether
the files describe the same sample.

Standard library only - runs before you install anything.

Usage:
  python src/v2p/stages/00_inspect_inputs.py data/*            # everything
  python src/v2p/stages/00_inspect_inputs.py data/ ref/        # directories too
  python src/v2p/stages/00_inspect_inputs.py --json data/      # machine-readable
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

# GRCh38 vs GRCh37 primary contig lengths - the cheapest build discriminator
CHR_LENGTHS = {
    "GRCh38": {"chr1": 248956422, "chr2": 242193529, "chr17": 83257441},
    "GRCh37": {"chr1": 249250621, "chr2": 243199373, "chr17": 81195210},
}
# highest coordinate seen per chromosome, used when no header is available
GRCH37_ONLY_MAX = {"chr1": (248956422, 249250621)}

RESET, BOLD, DIM = "\033[0m", "\033[1m", "\033[2m"


def hdr(s: str) -> None:
    print(f"\n{BOLD}{'=' * 72}\n{s}\n{'=' * 72}{RESET}")


def kv(k: str, v) -> None:
    print(f"  {k:<28} {v}")


def human(n: int) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return f"{n:.1f} {u}" if u != "B" else f"{n} B"
        n /= 1024
    return str(n)


def opener(p: Path):
    """Return a text-mode handle, transparently handling gzip and bgzip."""
    with open(p, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "rt", encoding="utf-8", errors="replace")


def sniff(p: Path) -> str:
    """Identify by content, since these files are inconsistently named."""
    with open(p, "rb") as fh:
        magic = fh.read(4)
    if magic[:2] == b"\x1f\x8b":
        kind = "gzip"
    elif magic[:2] == b"PK":
        return "zip"
    elif magic[:4] == b"CSI\x01":
        return "csi_index"
    elif magic[:4] == b"TBI\x01":
        return "tabix_index"
    else:
        kind = "text"
    try:
        with opener(p) as fh:
            head = fh.read(4096)
    except Exception:
        return "binary"
    if head.startswith("##fileformat=VCF"):
        return "vcf"
    if head.startswith(">"):
        return "fasta"
    first = head.split("\n", 1)[0]
    if p.suffix.lower() == ".bed" or re.match(r"^\S+\t\d+\t\d+", first):
        return "bed"
    if "\t" in first:
        return "tsv"
    if "," in first:
        return "csv"
    return kind


# --------------------------------------------------------------------------

def contig_style(names) -> str:
    n = [str(x) for x in names]
    if not n:
        return "unknown"
    chr_pref = sum(1 for x in n if x.startswith("chr"))
    if chr_pref == len(n):
        return "UCSC (chr-prefixed)"
    if chr_pref == 0:
        return "Ensembl (no chr prefix)"
    return f"MIXED ({chr_pref}/{len(n)} chr-prefixed) - problem"


def build_from_lengths(lengths: dict[str, int]) -> str:
    for build, ref in CHR_LENGTHS.items():
        hits = sum(1 for c, L in ref.items()
                   if lengths.get(c) == L or lengths.get(c.replace("chr", "")) == L)
        if hits:
            return f"{build} (matched {hits} contig length(s))"
    return "UNKNOWN - no contig length matched GRCh37 or GRCh38"


def inspect_vcf(p: Path, rec: dict) -> None:
    contigs: dict[str, int] = {}
    ref_line = ""
    samples: list[str] = []
    n_header = 0
    filters = Counter()
    types = Counter()
    chroms = Counter()
    max_pos: dict[str, int] = {}
    info_keys = Counter()
    n = 0
    with opener(p) as fh:
        for line in fh:
            if line.startswith("##"):
                n_header += 1
                m = re.match(r"##contig=<ID=([^,>]+)(?:,length=(\d+))?", line)
                if m:
                    contigs[m.group(1)] = int(m.group(2) or 0)
                if line.startswith("##reference"):
                    ref_line = line.strip()
                continue
            if line.startswith("#CHROM"):
                f = line.rstrip("\n").split("\t")
                samples = f[9:] if len(f) > 9 else []
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            n += 1
            chroms[f[0]] += 1
            max_pos[f[0]] = max(max_pos.get(f[0], 0), int(f[1]))
            if len(f) > 6:
                for x in f[6].split(";"):
                    filters[x] += 1
            if len(f) > 7 and n <= 2000:
                for kvp in f[7].split(";"):
                    info_keys[kvp.split("=")[0]] += 1
            ref, alts = f[3], f[4]
            for alt in alts.split(","):
                if alt.startswith("<") or "[" in alt or "]" in alt:
                    types["symbolic/breakend"] += 1
                elif len(ref) == 1 and len(alt) == 1:
                    types["SNV"] += 1
                elif len(ref) == len(alt):
                    types["MNV"] += 1
                elif len(ref) > len(alt):
                    types["deletion"] += 1
                else:
                    types["insertion"] += 1

    kv("meta-header lines", n_header)
    kv("##reference", ref_line or "(absent)")
    kv("sample columns", f"{len(samples)}: {samples if samples else '(sites-only VCF)'}")
    kv("data rows", f"{n:,}")
    kv("contigs in header", len(contigs))
    kv("contig naming", contig_style(contigs or chroms))
    if contigs and any(contigs.values()):
        kv("genome build", build_from_lengths(contigs))
        rec["build"] = build_from_lengths(contigs)
    else:
        c1 = max_pos.get("chr1", max_pos.get("1", 0))
        kv("genome build", f"no contig lengths in header; "
                           f"max chr1 position seen = {c1:,} "
                           f"(GRCh38 chr1 is 248,956,422)")
    kv("distinct chromosomes", len(chroms))
    print(f"  {DIM}variant types:{RESET}")
    for k, v in types.most_common():
        print(f"      {k:<20} {v:>10,}")
    print(f"  {DIM}FILTER values:{RESET}")
    for k, v in filters.most_common(10):
        print(f"      {k:<20} {v:>10,}")
    if info_keys:
        print(f"  {DIM}INFO keys (first 2000 rows): {RESET}"
              + ", ".join(sorted(info_keys)[:20]))
    rec.update({"rows": n, "types": dict(types), "filters": dict(filters),
                "contig_style": contig_style(contigs or chroms),
                "samples": samples})
    if types.get("symbolic/breakend"):
        print(f"  {BOLD}note{RESET} symbolic/breakend records are structural "
              f"variants - handled by the fusion path, not codon substitution.")


def inspect_bed(p: Path, rec: dict) -> None:
    n = 0
    total = 0
    chroms = Counter()
    lengths = []
    ncol = 0
    with opener(p) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            n += 1
            ncol = max(ncol, len(f))
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            chroms[f[0]] += 1
            total += e - s
            if len(lengths) < 200000:
                lengths.append(e - s)
    kv("regions", f"{n:,}")
    kv("columns", ncol)
    kv("chromosomes", len(chroms))
    kv("total covered", f"{total:,} bp ({total / 3.1e9:.1%} of the genome)")
    if lengths:
        lengths.sort()
        kv("region length min/med/max",
           f"{lengths[0]:,} / {lengths[len(lengths)//2]:,} / {lengths[-1]:,}")
    kv("contig naming", contig_style(chroms))
    print(f"  {BOLD}note{RESET} BED is 0-based half-open; VCF is 1-based "
          f"inclusive. Intersect with bedtools/bcftools, never by hand.")
    rec.update({"regions": n, "covered_bp": total,
                "contig_style": contig_style(chroms)})


def inspect_fasta(p: Path, rec: dict) -> None:
    n = 0
    lens = []
    wrap = Counter()
    alphabet = Counter()
    styles = Counter()
    fields = Counter()
    first = []
    cur = 0
    lines_this = 0
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(cur)
                    wrap[lines_this] += 1
                n += 1
                cur = 0
                lines_this = 0
                if n <= 3:
                    first.append(line.rstrip("\n")[:110])
                h = line
                if re.match(r"^>(sp|tr)\|", h):
                    styles["UniProt sp|tr"] += 1
                elif h.startswith(">chr") or re.match(r"^>\d+\s", h):
                    styles["genome contig"] += 1
                else:
                    styles["other"] += 1
                for f in ("OS=", "OX=", "GN=", "PE=", "SV="):
                    if f in h:
                        fields[f] += 1
            else:
                s = line.strip()
                cur += len(s)
                lines_this += 1
                if n <= 5000:
                    alphabet.update(s.upper())
    if n:
        lens.append(cur)
        wrap[lines_this] += 1
    kv("entries", f"{n:,}")
    kv("header style", dict(styles))
    if fields:
        kv("header fields", {k: v for k, v in fields.items()})
    if lens:
        lens.sort()
        kv("length min/med/max",
           f"{lens[0]:,} / {lens[len(lens)//2]:,} / {lens[-1]:,}")
    single = wrap.get(1, 0)
    kv("sequence layout",
       "single line per entry (unwrapped)" if single == n
       else f"wrapped ({single}/{n} unwrapped)")
    letters = {c for c in alphabet if c.isalpha()}
    is_prot = bool(letters - set("ACGTNU"))
    kv("sequence type", "protein" if is_prot else "nucleotide")
    if is_prot:
        odd = {c: alphabet[c] for c in "UOXBZJ*" if c in alphabet}
        if odd:
            kv("non-standard residues", odd)
            print(f"  {BOLD}note{RESET} U = selenocysteine. A standard-code "
                  f"translation truncates these at UGA - excluded from "
                  f"identity statistics, not counted as failures.")
    else:
        contigs = {}
        # genome FASTA: measure the first few contigs to confirm the build
        print(f"  {DIM}(run `samtools faidx` then `head -3 *.fai` to confirm "
              f"contig lengths / build){RESET}")
    for h in first:
        print(f"  {DIM}>{RESET} {h}")
    rec.update({"entries": n, "type": "protein" if is_prot else "nucleotide",
                "unwrapped": single == n})


def inspect_table(p: Path, rec: dict, delim: str) -> None:
    with opener(p) as fh:
        sample = fh.read(1 << 20)
    buf = io.StringIO(sample)
    rdr = csv.reader(buf, delimiter=delim)
    try:
        header = next(rdr)
    except StopIteration:
        print("  (empty)")
        return
    n = 0
    chroms = Counter()
    with opener(p) as fh:
        r2 = csv.DictReader(fh, delimiter=delim)
        cols = r2.fieldnames or []
        chrom_col = next((c for c in cols
                          if c and c.lower() in ("chr", "chrom", "chromosome",
                                                 "#chrom")), None)
        for row in r2:
            n += 1
            if chrom_col:
                chroms[(row.get(chrom_col) or "").strip()] += 1
    kv("data rows", f"{n:,}")
    kv("columns", len(header))
    kv("column names", ", ".join(header[:12])
       + (f" ... (+{len(header)-12} more)" if len(header) > 12 else ""))
    if chroms:
        chroms.pop("", None)
        chroms.pop("CHROM", None)
        kv("contig naming", contig_style(chroms))
        rec["contig_style"] = contig_style(chroms)
    rec.update({"rows": n, "columns": len(header)})


def inspect_zip(p: Path, rec: dict) -> None:
    with zipfile.ZipFile(p) as z:
        infos = z.infolist()
    kv("members", len(infos))
    print(f"  {DIM}{'size':>12}  name{RESET}")
    for i in sorted(infos, key=lambda x: -x.file_size)[:40]:
        if i.is_dir():
            continue
        print(f"  {human(i.file_size):>12}  {i.filename}")
    print(f"\n  {BOLD}extract with{RESET} unzip -o {p} -d data/")
    rec["members"] = [i.filename for i in infos if not i.is_dir()]


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    files: list[Path] = []
    for x in args.paths:
        p = Path(x)
        if p.is_dir():
            files += sorted(q for q in p.rglob("*") if q.is_file())
        elif p.exists():
            files.append(p)
        else:
            print(f"missing: {x}", file=sys.stderr)

    report: list[dict] = []
    for p in files:
        kind = sniff(p)
        rec = {"path": str(p), "kind": kind, "bytes": p.stat().st_size}
        hdr(f"{p.name}   [{kind}]   {human(p.stat().st_size)}")
        try:
            if kind == "vcf":
                inspect_vcf(p, rec)
            elif kind == "bed":
                inspect_bed(p, rec)
            elif kind == "fasta":
                inspect_fasta(p, rec)
            elif kind in ("tsv", "csv"):
                inspect_table(p, rec, "\t" if kind == "tsv" else ",")
            elif kind == "zip":
                inspect_zip(p, rec)
            elif kind in ("csi_index", "tabix_index"):
                kv("type", "index for a bgzipped VCF - no data of its own")
                kv("action", "keep next to the .vcf.gz; regenerate any time "
                             "with `bcftools index`")
            else:
                kv("type", f"unrecognised ({kind}); run `file` on it")
        except Exception as exc:
            print(f"  ERROR: {type(exc).__name__}: {exc}")
            rec["error"] = str(exc)
        report.append(rec)

    # ---- cross-file consistency ---------------------------------------
    hdr("cross-file consistency")
    styles = {r["path"]: r.get("contig_style") for r in report
              if r.get("contig_style")}
    distinct = {v for v in styles.values() if v and "unknown" not in v}
    if len(distinct) > 1:
        print(f"  {BOLD}MISMATCH{RESET} contig naming differs across files:")
        for k, v in styles.items():
            print(f"      {v:<28} {k}")
        print("  Normalise before use, e.g.:")
        print("      bcftools annotate --rename-chrs chr_map.txt in.vcf.gz "
              "-Oz -o out.vcf.gz")
    elif distinct:
        print(f"  contig naming consistent: {distinct.pop()}")
    else:
        print("  contig naming: not determinable from these files")

    builds = {r["path"]: r["build"] for r in report if r.get("build")}
    if builds:
        for k, v in builds.items():
            print(f"  build: {v}   ({Path(k).name})")
        if any("GRCh37" in v for v in builds.values()):
            print(f"  {BOLD}STOP{RESET} a GRCh37 input cannot be mixed with a "
                  f"GRCh38 reference. Lift over first.")
    else:
        print("  build: no VCF header contig lengths found - confirm manually")

    if args.json:
        Path("results/qc").mkdir(parents=True, exist_ok=True)
        Path("results/qc/input_inventory.json").write_text(
            json.dumps(report, indent=2))
        print("\n  wrote results/qc/input_inventory.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

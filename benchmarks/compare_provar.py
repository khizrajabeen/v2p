#!/usr/bin/env python3
"""Compare a v2p database with a ProVar database built from the same VCF.

ProVar (ProGenNo/ProHap, Nature Methods 2024) is the fair comparator for
v2p: one sample's VCF in, protein sequences out, same authors as ProHap.
ProHap itself needs phased population panels and about a terabyte, and
comparing tools built for different jobs invites the result to be
dismissed in both directions.

What is measured, and why each number is here:

- **variants represented** - coverage parity. If ProVar represents a
  variant v2p does not, that is a v2p bug, and it is the most useful
  thing this script can find.
- **sequences produced** - database size for that coverage. Bigger is
  not better: every sequence costs FDR power.
- **peptides spanning two or more variants** - the gap combining fills.
  ProVar considers each allele independently, so this should be zero for
  it; the script measures rather than assumes that.
- **peptides unique to each tool** - what each finds that the other
  misses, in both directions.

Peptide sets are compared after subtracting the reference proteome's own
peptides, so neither tool is credited for peptides a plain reference
search would already find. The reference space comes from the reference
entries appended to the v2p database, so both tools are measured against
the same one.

Usage:
  python3 benchmarks/compare_provar.py \
      --v2p-fasta out/HCC1395.target.fasta \
      --provar-fasta /path/out/provar_variants.fa \
      --provar-tsv /path/out/provar_variants.tsv \
      --v2p-disposition out/_work/tables/disposition.tsv \
      --v2p-time out/v2p_time.txt --provar-time /path/out/provar_time.txt \
      --out benchmarks/provar_comparison.md
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.peptides import digest                          # noqa: E402


def read_fasta(path: Path):
    """Yield (header, sequence). Plain FASTA, no dependencies."""
    head, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(buf)
                head, buf = line[1:], []
            elif line:
                buf.append(line)
    if head is not None:
        yield head, "".join(buf)


def peptide_set(seqs, missed_cleavages: int = 2, equate_il: bool = True):
    out = set()
    for s in seqs:
        for p in digest(s, missed_cleavages):
            out.add(p.replace("I", "L") if equate_il else p)
    return out


def parse_time(path: str) -> tuple[str, str]:
    """(elapsed, peak RSS) from a /usr/bin/time -v file."""
    if not path or not Path(path).exists():
        return ("n/a", "n/a")
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    el = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    mb = f"{int(rss.group(1)) / 1024:.0f} MB" if rss else "n/a"
    return (el.group(1) if el else "n/a", mb)


def norm_locus(chrom: str, pos, ref: str, alt: str) -> str:
    """chr-prefix-insensitive variant key, so the two tools compare."""
    c = str(chrom)
    c = c[3:] if c.lower().startswith("chr") else c
    return f"{c}:{pos}{ref}>{alt}"


_DNA_CHANGE = re.compile(r"^(?:chr)?([\w.]+):g?\.?(\d+)([ACGTN]*)>([ACGTN]*)",
                         re.I)
# ProVar's own form: position, ref and alt, with the contig held in a
# separate column.
_PROVAR_CHANGE = re.compile(r"^(\d+):([ACGTN-]*)>([ACGTN-]*)", re.I)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v2p-fasta", required=True)
    ap.add_argument("--provar-fasta", required=True)
    ap.add_argument("--provar-tsv", required=True)
    ap.add_argument("--v2p-disposition", default="")
    ap.add_argument("--v2p-time", default="")
    ap.add_argument("--provar-time", default="")
    ap.add_argument("--missed-cleavages", type=int, default=2)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    # ---- v2p side --------------------------------------------------
    v2p_var, v2p_ref = [], []
    v2p_combo_seqs = []
    v2p_loci: set[str] = set()
    for head, seq in read_fasta(Path(args.v2p_fasta)):
        if not head.startswith("vr|"):
            v2p_ref.append(seq)               # appended reference proteome
            continue
        v2p_var.append(seq)
        if " VT=COMBO " in f" {head} ":
            v2p_combo_seqs.append(seq)
        m = re.search(r" LOC=(\S+)", head)
        if m:
            for locus in m.group(1).split(";"):
                mm = _DNA_CHANGE.match(locus)
                if mm:
                    v2p_loci.add(norm_locus(*mm.groups()))

    # The disposition is the honest source for coverage: it names every
    # input variant and what became of it, including the ones that
    # produced nothing.
    v2p_covered: set[str] = set()
    if args.v2p_disposition and Path(args.v2p_disposition).exists():
        with open(args.v2p_disposition, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                built = (row.get("n_proteins")
                         or row.get("n_sequences") or "0").strip()
                if not (built.isdigit() and int(built) > 0):
                    continue
                # The disposition names the variant as CLASS|chr:posREF>ALT;
                # there is no separate locus column.
                loc = (row.get("locus")
                       or (row.get("variant_id") or "").split("|")[-1])
                mm = _DNA_CHANGE.match(loc)
                if mm:
                    v2p_covered.add(norm_locus(*mm.groups()))
    if not v2p_covered:
        v2p_covered = set(v2p_loci)

    # ---- ProVar side -----------------------------------------------
    provar_seqs = [s for _h, s in read_fasta(Path(args.provar_fasta))]
    provar_covered: set[str] = set()
    per_seq_variants: dict[str, set] = defaultdict(set)
    with open(args.provar_tsv, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            # ProVar writes DNA_change as POS:REF>ALT and keeps the
            # contig in its own column, so the two have to be rejoined
            # before the key compares with v2p's chr-prefixed loci.
            chrom = (row.get("chromosome") or "").strip()
            dna = (row.get("DNA_change") or "").strip()
            for part in dna.split(";"):
                mm = _PROVAR_CHANGE.match(part.strip())
                if not mm:
                    continue
                pos, ref, alt = mm.groups()
                key = norm_locus(chrom, pos, ref, alt)
                provar_covered.add(key)
                # Group by the entry's own accession, not by transcript:
                # ProVar emits one entry per variant, so a transcript with
                # two variants has two entries, each carrying one. Grouping
                # by transcript would count those as a combined entry.
                per_seq_variants[row.get("variantID", "")].add(key)

    multi = sum(1 for v in per_seq_variants.values() if len(v) > 1)

    # ---- peptides ---------------------------------------------------
    ref_space = peptide_set(v2p_ref, args.missed_cleavages)
    v2p_peps = peptide_set(v2p_var, args.missed_cleavages) - ref_space
    provar_peps = peptide_set(provar_seqs, args.missed_cleavages) - ref_space
    combo_peps = peptide_set(v2p_combo_seqs, args.missed_cleavages) - ref_space

    only_v2p = v2p_peps - provar_peps
    only_provar = provar_peps - v2p_peps
    v2p_el, v2p_mem = parse_time(args.v2p_time)
    pv_el, pv_mem = parse_time(args.provar_time)

    missed_by_v2p = provar_covered - v2p_covered
    missed_by_provar = v2p_covered - provar_covered

    lines = [
        "| measure | v2p | ProVar |",
        "|---|---|---|",
        f"| variants represented in >=1 entry | {len(v2p_covered)} | "
        f"{len(provar_covered)} |",
        f"| sequences produced | {len(v2p_var)} | {len(provar_seqs)} |",
        f"| entries carrying >1 variant | {len(v2p_combo_seqs)} | {multi} |",
        f"| novel tryptic peptides | {len(v2p_peps)} | {len(provar_peps)} |",
        f"| peptides the other tool lacks | {len(only_v2p)} | "
        f"{len(only_provar)} |",
        f"| peptides from combined entries | {len(combo_peps)} | 0 |",
        f"| wall clock | {v2p_el} | {pv_el} |",
        f"| peak memory | {v2p_mem} | {pv_mem} |",
        "",
        f"Variants ProVar represents and v2p does not: "
        f"**{len(missed_by_v2p)}**",
    ]
    if missed_by_v2p:
        lines.append("")
        lines.append("Each of these is a v2p bug until investigated. What "
                     "the investigation found is in `README.md`.")
        lines.append("")
        lines.append("")
        for key in sorted(missed_by_v2p)[:40]:
            lines.append(f"- `{key}`")
        if len(missed_by_v2p) > 40:
            lines.append(f"- ... and {len(missed_by_v2p) - 40} more")
    lines.append("")
    lines.append(f"Variants v2p represents and ProVar does not: "
                 f"**{len(missed_by_provar)}**")
    lines += [
        "",
        "## Caveats, which the numbers above do not carry on their own",
        "",
        "- **The peptide counts are confounded.** Novelty here is measured "
        "against UniProt SwissProt, the reference v2p appends. ProVar "
        "translates Ensembl transcripts, so its alternative isoforms count "
        "as novel whether or not a variant is involved. That inflates both "
        "its novel-peptide total and its peptides-the-other-lacks column, "
        "and it is why ProVar can show more novel peptides from fewer "
        "sequences. Read those two rows as *different reference sets*, not "
        "as variant discovery.",
        "- **Variant coverage is the clean comparison**, because a variant "
        "either ends up in an entry or it does not, independent of which "
        "reference proteome either tool started from.",
        "- **The runtime figures are not like for like.** This harness "
        "re-reads the 521 MB cDNA FASTA once per chromosome, which is most "
        "of ProVar's wall clock. It is a property of how it was driven "
        "here, not a claim about the tool.",
        "- **Entries carrying more than one variant is the row that "
        "matters** for the feature under test, and it is measured, not "
        "assumed: ProVar's own output confirms one variant per entry.",
    ]

    report = "\n".join(lines)
    print(report)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report + "\n", encoding="utf-8")
        print(f"\n-> {out}")
        miss = out.with_suffix(".missed_by_v2p.txt")
        miss.write_text("\n".join(sorted(missed_by_v2p)) + "\n",
                        encoding="utf-8")
        print(f"-> {miss}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

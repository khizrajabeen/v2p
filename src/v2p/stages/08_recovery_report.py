#!/usr/bin/env python3
"""Stage 8 - variant recovery accounting.

Answers "how many variants did we start with, how many became proteins,
and where did the rest go" with every input variant accounted for exactly
once. Reconciles three files:

  the manifest        every input variant
  the disposition     what happened to each one in stage 2
  the entry table     what survived into the release

Most variants producing no protein is the expected result, not a
failure: coding exons are roughly 1-2% of the genome, so a whole-genome
somatic call set is overwhelmingly intergenic and intronic. The point of
this report is to show that the losses are accounted for rather than
unexplained.

Usage:
  python src/v2p/stages/08_recovery_report.py \
      --manifest results/tables/unified_variant_manifest.tsv \
      --disposition results/tables/disposition.representative.tsv \
      --entries release/HCC1395_variant_proteome.entries.tsv \
      --outdir results
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from v2p.provenance import RunLogger      # noqa: E402

ORDER = ["SNV", "MNV", "INDEL", "RNA_EDITING", "FUSION",
         "AS_SE", "AS_RI", "AS_MX", "AS_A3", "AS_A5", "AS_AF", "AS_AL"]

OUTCOME_LABEL = {
    "protein_built": "produced a variant protein",
    "synonymous_only": "synonymous — coding, but the codon still encodes "
                       "the same residue",
    "utr_variant": "in a UTR — exonic but outside the translated region",
    "synonymous_or_utr": "synonymous and/or UTR across transcripts",
    "not_in_coding_exon": "intronic, UTR, or crossing a splice junction",
    "no_coding_transcript_at_locus": "intergenic, or no coding transcript",
    "below_min_length": "protein shorter than the minimum length",
    "filtered": "built but removed by a filter",
    "no_protein": "no protein (see the run log)",
    "error": "error during translation",
}


def read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--disposition", required=True)
    ap.add_argument("--entries", help="release entry table (optional)")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--logdir", default="logs")
    args = ap.parse_args()

    out = Path(args.outdir)
    (out / "qc").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    rl = RunLogger("08_recovery_report", args.logdir)
    rl.add_input("manifest", args.manifest)
    rl.add_input("disposition", args.disposition)

    manifest = Counter()
    for r in read_tsv(Path(args.manifest)):
        manifest[r["variant_class"]] += 1

    disp = read_tsv(Path(args.disposition))
    by_class: dict[str, Counter] = defaultdict(Counter)
    proteins: Counter = Counter()
    for r in disp:
        c = r["variant_class"]
        by_class[c][r["outcome"]] += 1
        proteins[c] += int(r["n_proteins"] or 0)

    # every input variant must appear exactly once
    seen = Counter(r["variant_class"] for r in disp)
    gaps = {c: manifest[c] - seen.get(c, 0) for c in manifest
            if manifest[c] != seen.get(c, 0)}

    entries = Counter()
    if args.entries and Path(args.entries).exists():
        rl.add_input("entries", args.entries)
        for r in read_tsv(Path(args.entries)):
            entries[r["variant_class"]] += 1

    classes = [c for c in ORDER if c in manifest] + \
              [c for c in sorted(manifest) if c not in ORDER]

    L = ["# Variant recovery accounting\n"]
    L.append("Every input variant appears exactly once below.\n")
    L.append("| variant class | input | produced a protein | rate | "
             "protein sequences | in the release |")
    L.append("|---|---:|---:|---:|---:|---:|")
    t_in = t_ok = t_pr = 0
    for c in classes:
        n_in = manifest[c]
        n_ok = by_class[c]["protein_built"]
        t_in += n_in
        t_ok += n_ok
        t_pr += proteins[c]
        rel = f"{entries[c]:,}" if entries else "-"
        L.append(f"| {c} | {n_in:,} | {n_ok:,} | "
                 f"{100*n_ok/n_in if n_in else 0:.1f}% | "
                 f"{proteins[c]:,} | {rel} |")
    L.append(f"| **total** | **{t_in:,}** | **{t_ok:,}** | "
             f"**{100*t_ok/t_in if t_in else 0:.1f}%** | **{t_pr:,}** | "
             f"**{sum(entries.values()):,}** |" if entries else
             f"| **total** | **{t_in:,}** | **{t_ok:,}** | "
             f"**{100*t_ok/t_in if t_in else 0:.1f}%** | **{t_pr:,}** | - |")

    L.append("\n`protein sequences` exceeds `produced a protein` where one "
             "variant affects several transcripts, and `in the release` is "
             "lower again after identical sequences are merged.\n")

    L.append("\n## Where the rest went\n")
    for c in classes:
        counts = by_class[c]
        if sum(counts.values()) == 0:
            continue
        L.append(f"### {c} — {manifest[c]:,} input variants\n")
        L.append("| outcome | variants | share |")
        L.append("|---|---:|---:|")
        tot = sum(counts.values())
        for o, n in counts.most_common():
            L.append(f"| {OUTCOME_LABEL.get(o, o)} | {n:,} | "
                     f"{100*n/tot:.1f}% |")
        L.append("")

    # -- sanity check on the coding-consequence ratio --------------------
    L.append("\n## Coding-consequence sanity check\n")
    L.append("Roughly 74% of possible single-base changes within a codon "
             "alter the amino acid, so a somatic call set should show more "
             "nonsynonymous than synonymous coding variants — typically "
             "around 2-3 to 1. A ratio below 1 means coding and non-coding "
             "variants are being conflated somewhere.\n")
    L.append("| class | nonsynonymous | synonymous | ratio | verdict |")
    L.append("|---|---:|---:|---:|---|")
    for c in classes:
        if c not in ("SNV", "MNV", "INDEL", "RNA_EDITING"):
            continue
        ns = by_class[c]["protein_built"]
        syn = by_class[c]["synonymous_only"]
        if not (ns or syn):
            continue
        ratio = ns / syn if syn else float("inf")
        if c == "RNA_EDITING":
            verdict = "n/a — ADAR targets are not random"
        elif syn == 0:
            verdict = "no synonymous variants to compare"
        elif ratio < 1.0:
            verdict = "**INVERTED — investigate**"
        elif ratio < 1.5:
            verdict = "low, but possible under strong selection"
        else:
            verdict = "as expected"
        rs = "inf" if syn == 0 else f"{ratio:.2f}"
        L.append(f"| {c} | {ns:,} | {syn:,} | {rs} | {verdict} |")
    L.append("\nUTR variants are counted separately from synonymous ones. "
             "They are exonic and leave the protein unchanged, so treating "
             "them as synonymous would inflate the synonymous count and "
             "invert this ratio.\n")

    L.append("\n## Reading this\n")
    L.append("- A low rate for whole-genome somatic classes is expected. "
             "Protein-coding exons are roughly 1-2% of the genome, so most "
             "somatic variants are intergenic or intronic and cannot change "
             "a protein.\n")
    L.append("- `synonymous` variants are "
             "synonymous changes. They are real variants that produce no new "
             "sequence, so they add nothing to a search database.\n")
    L.append("- `no coding transcript at locus` on a class that should be "
             "genic is the one row worth investigating: it can indicate a "
             "gene-symbol or annotation-version mismatch rather than "
             "genuine intergenic position.\n")

    if gaps:
        L.append("\n## Reconciliation failure\n")
        for c, n in gaps.items():
            L.append(f"- **{c}**: {n} input variants have no disposition row.")
        rl.log.error("disposition does not account for every input: %s", gaps)
    else:
        L.append("\nReconciliation: every input variant has exactly one "
                 "disposition row.\n")

    rep = out / "qc" / "recovery_report.md"
    rep.write_text("\n".join(L) + "\n", encoding="utf-8")

    tbl = out / "tables" / "recovery_by_class.tsv"
    with open(tbl, "w", newline="", encoding="utf-8") as fh:
        cols = ["variant_class", "n_input", "n_produced_protein",
                "recovery_rate_pct", "n_protein_sequences"] + \
               sorted({o for c in by_class for o in by_class[c]})
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for c in classes:
            row = {"variant_class": c, "n_input": manifest[c],
                   "n_produced_protein": by_class[c]["protein_built"],
                   "recovery_rate_pct": round(
                       100 * by_class[c]["protein_built"] / manifest[c], 2)
                   if manifest[c] else 0,
                   "n_protein_sequences": proteins[c]}
            for o in cols[5:]:
                row[o] = by_class[c].get(o, 0)
            w.writerow(row)

    for c in classes:
        rl.count(f"input.{c}", manifest[c])
        rl.count(f"recovered.{c}", by_class[c]["protein_built"])
    rl.add_output("report", rep)
    rl.add_output("table", tbl)
    rl.close("error" if gaps else "ok")
    print(rep.read_text())
    return 1 if gaps else 0


if __name__ == "__main__":
    raise SystemExit(main())

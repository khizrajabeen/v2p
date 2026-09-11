#!/usr/bin/env python3
"""A worked cross-evidence proteoform: an ADAR edit plus a somatic SNV.

**This is a constructed case, not a measurement.** Read that first.

The claim it illustrates is the one thing v2p does that a genotype-based
tool cannot: combining an A-to-I RNA editing site with a somatic DNA
variant on one transcript. No haplotype panel can represent that
proteoform, because RNA editing is not in the genome and never appears
in a VCF of genotypes.

HCC1395 does not contain an example. That was checked, not assumed, at
three levels:

  * transcripts carrying both an editing site and a DNA variant: **40**
  * of those, transcripts where *both* recode: **0** (the edits there
    are synonymous or fall in UTR)
  * genes with a recoding edit AND a recoding DNA variant, even on
    different transcripts: **0** (29 genes have the first, 278 the
    second, and the sets do not intersect)

So one half of this case is real and one half is invented, and which is
which matters:

  * the **RNA editing site is real** - a called A-to-I site from
    HCC1395's own RES table, on a transcript where it genuinely recodes.
    This is the half no other tool can use.
  * the **SNV is synthetic** - placed in the same tryptic peptide so the
    two share a peptide. It stands in for a somatic call that HCC1395
    happens not to have at that position.

What this demonstrates is the mechanism and the code path. What it does
not demonstrate is frequency in nature; for that, a sample with both
call types co-located is needed, and the search for one is recorded in
benchmarks/README.md.

Usage:
  python3 benchmarks/cross_evidence_case.py \
      --manifest results/tables/unified_variant_manifest.tsv \
      --genome ref/GRCh38.primary_assembly.genome.fa \
      --gtf ref/gencode.v44.annotation.gtf.gz \
      --out benchmarks/cross_evidence_case.md
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.annotation import Annotation, Genome              # noqa: E402
from v2p.build.combinatorial import (                      # noqa: E402
    CombinableVariant, _placement, _translate,
    build_combinatorial_proteins, cooccurring_peptides,
)
from v2p.build.smallvar import _cds_offset, _tx_sequence   # noqa: E402
from v2p.seqops import table_for_contig, translate_orf     # noqa: E402

SWAP = {"A": "C", "C": "A", "G": "T", "T": "G"}


def peptide_bounds(prot: str, idx: int, missed: int = 2):
    """(start, end) of the longest tryptic peptide containing `idx`."""
    cuts = [0] + [m.start() for m in
                  re.finditer(r"(?<=[KR])(?!P)", prot)] + [len(prot)]
    best = None
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + missed + 2, len(cuts))):
            s, e = cuts[i], cuts[j]
            if s <= idx < e and 6 <= e - s <= 50:
                if best is None or (e - s) > (best[1] - best[0]):
                    best = (s, e)
    return best


def genomic_of_codon(t, cds: int, codon_idx: int):
    """Genomic position of the first base of a codon on transcript `t`."""
    tx_off = cds + 3 * codon_idx
    seen = 0
    blocks = t.exons if t.strand == "+" else list(reversed(t.exons))
    for s, e in blocks:
        ln = e - s + 1
        if tx_off < seen + ln:
            return (s + (tx_off - seen) if t.strand == "+"
                    else e - (tx_off - seen))
        seen += ln
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    ann = Annotation.from_gtf(Path(args.gtf))
    genome = Genome(Path(args.genome))

    edits = []
    with open(args.manifest, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["variant_class"] != "RNA_EDITING":
                continue
            p = json.loads(row["payload_json"])
            if len(p["ref"]) == 1 and len(p["alt"]) == 1:
                edits.append(CombinableVariant(
                    p["chrom"], p["pos"], p["ref"], p["alt"],
                    "RNA_EDITING", row["source"]))

    chosen = None
    for ed in edits:
        for t in ann.representatives_at(ed.chrom, ed.pos):
            if not t.is_coding:
                continue
            place = _placement(t, ed)
            cds = _cds_offset(t)
            if place is None or cds is None or place[0] < cds:
                continue
            tx_seq = _tx_sequence(t, genome)
            ref_p = translate_orf(
                tx_seq[cds:], table_for_contig(t.chrom),
                sec_codons=t.sec_codon_indices(cds)).protein
            alt_p = _translate(tx_seq, [place], cds, t)
            if not alt_p or alt_p == ref_p:
                continue                      # synonymous: not recoding
            idx = (place[0] - cds) // 3
            bounds = peptide_bounds(ref_p, idx)
            if bounds is None:
                continue
            chosen = (ed, t, cds, tx_seq, ref_p, alt_p, idx, bounds)
            break
        if chosen:
            break

    if not chosen:
        print("no recoding editing site with a usable peptide found")
        return 1

    ed, t, cds, tx_seq, ref_p, alt_p, idx, (pstart, pend) = chosen

    # A synthetic SNV in the same peptide, at a codon that is not the
    # edited one. Its genomic position is derived from the transcript
    # offset, so the coordinates are real even though the call is not.
    target = None
    for cand in range(pstart, pend):
        if cand == idx:
            continue
        gpos = genomic_of_codon(t, cds, cand)
        if gpos is None:
            continue
        base = genome.fetch(t.chrom, gpos, gpos).upper()
        if base not in SWAP:
            continue
        snv = CombinableVariant(t.chrom, gpos, base, SWAP[base], "SNV",
                                "synthetic")
        pl = _placement(t, snv)
        if pl is None:
            continue
        only_snv = _translate(tx_seq, [pl], cds, t)
        if only_snv and only_snv != ref_p:
            target = (snv, only_snv, cand)
            break

    if not target:
        print("could not place a synthetic SNV in the same peptide")
        return 1
    snv, only_snv, snv_idx = target

    recs = build_combinatorial_proteins([ed, snv], ann, genome,
                                        transcript_mode="representative")
    if not recs:
        print("the combination produced no entry")
        return 1
    rec = recs[0]
    peps = cooccurring_peptides(rec.sequence, [alt_p, only_snv], ref_p)
    spanning = sorted(peps, key=len)[:3]

    gene = t.gene_name or t.gene_id
    lo, hi = min(idx, snv_idx), max(idx, snv_idx)
    lines = [
        "# A cross-evidence proteoform, constructed",
        "",
        "**This is a constructed case, not a measurement.** One half is a "
        "real call from HCC1395; the other is synthetic, and which is "
        "which is stated below.",
        "",
        "## Why it had to be constructed",
        "",
        "HCC1395 contains no transcript where an RNA editing site and a "
        "DNA variant both recode. That was measured, not assumed:",
        "",
        "| check | count |",
        "|---|---|",
        "| transcripts carrying both an editing site and a DNA variant | 40 |",
        "| of those, transcripts where **both** recode | **0** |",
        "| genes with a recoding edit *and* a recoding DNA variant | **0** |",
        "",
        "29 genes carry a recoding edit and 278 carry a recoding DNA "
        "variant; the two sets do not intersect.",
        "",
        "## The case",
        "",
        f"Gene **{gene}**, transcript `{t.tx_id}`.",
        "",
        "| element | status | detail |",
        "|---|---|---|",
        f"| A-to-I editing site | **real**, called in HCC1395 | "
        f"`{ed.locus}`, residue {idx + 1} |",
        f"| somatic SNV | **synthetic** | `{snv.locus}`, residue "
        f"{snv_idx + 1} |",
        "",
        "The editing site is the half no genotype-based tool can use. The "
        "SNV stands in for a somatic call this sample does not have at "
        "that position.",
        "",
        "## What the entries contain",
        "",
        f"| entry | residue {lo + 1} | residue {hi + 1} |",
        "|---|---|---|",
        f"| reference | {ref_p[lo]} | {ref_p[hi]} |",
        f"| edit only | {alt_p[lo]} | {alt_p[hi]} |",
        f"| SNV only | {only_snv[lo]} | {only_snv[hi]} |",
        f"| **combined** | {rec.sequence[lo]} | {rec.sequence[hi]} |",
        "",
        f"v2p describes the combined entry as `{rec.protein_change}`, "
        f"class `{rec.variant_class}`, notes `{','.join(rec.notes)}`.",
        "",
        "## The peptide that exists only in the combined form",
        "",
        "```",
    ]
    lines += list(spanning) or ["(none)"]
    lines += [
        "```",
        "",
        f"{len(peps)} tryptic peptide(s) of the combined protein appear in "
        "neither single-variant entry nor the reference. A search against "
        "a database built one variant at a time cannot identify them at "
        "any FDR.",
        "",
        "## What this does and does not show",
        "",
        "It shows the mechanism, and exercises the real code path on a "
        "real transcript with a real editing site. It does **not** show "
        "how often such proteoforms occur in nature. That needs a sample "
        "with both call types co-located, which the search recorded in "
        "`benchmarks/README.md` did not find.",
    ]

    report = "\n".join(lines)
    print(report)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report + "\n", encoding="utf-8")
        print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

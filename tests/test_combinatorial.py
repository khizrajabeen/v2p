#!/usr/bin/env python3
"""Combinatorial proteoforms - co-occurring variants applied together.

Run:  python tests/test_combinatorial.py
Exit code 0 = all pass. No pytest dependency so it runs anywhere.

The load-bearing assertion is that combining two variants yields a protein
carrying *both* changes, which neither single-variant entry contains. That
is the whole point: a peptide spanning both exists only in the combined
form, so a single-variant database cannot identify it at any FDR.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from make_fixture import build                              # noqa: E402
from v2p.annotation import Annotation, Genome               # noqa: E402
from v2p.build.combinatorial import (                       # noqa: E402
    COMBO_CLASS, CombinableVariant, build_combinatorial_proteins,
    group_by_transcript,
)
from v2p.build.smallvar import build_small_variant_proteins  # noqa: E402

FIX = ROOT / "tests" / "fixtures"
PASS: list[str] = []
FAIL: list[str] = []

SWAP = {"A": "C", "C": "A", "G": "T", "T": "G"}


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{' | ' + detail if detail else ''}")


def gpos_plus(meta, gene: str, codon: int, base: int = 1) -> int:
    m = meta[gene] if gene in meta else next(
        v for k, v in meta.items() if k.split("@")[0] == gene)
    tx_off = m["utr5"] + 3 * (codon - 1) + (base - 1)
    seen = 0
    for s, e in m["exons"]:
        ln = e - s + 1
        if tx_off < seen + ln:
            return s + (tx_off - seen)
        seen += ln
    raise IndexError


def main() -> int:
    meta = build()
    ann = Annotation.from_gtf(FIX / "mini.gtf")
    gen = Genome(FIX / "mini.fa")

    # Two substitutions on one transcript, several codons apart.
    p3 = gpos_plus(meta, "GPLUS", 3, 1)
    p9 = gpos_plus(meta, "GPLUS", 9, 1)
    r3 = gen.fetch("chrT1", p3, p3)
    r9 = gen.fetch("chrT1", p9, p9)
    a3, a9 = SWAP[r3], SWAP[r9]

    v3 = CombinableVariant("chrT1", p3, r3, a3, "SNV", "vcf")
    v9 = CombinableVariant("chrT1", p9, r9, a9, "SNV", "vcf")

    # What the single-variant path produces, for comparison.
    s3 = build_small_variant_proteins("chrT1", p3, r3, a3, ann, gen,
                                      variant_class="SNV", genes="GPLUS")
    s9 = build_small_variant_proteins("chrT1", p9, r9, a9, ann, gen,
                                      variant_class="SNV", genes="GPLUS")
    check("both single variants translate on their own",
          bool(s3) and bool(s9), f"{len(s3)}/{len(s9)}")

    # ------------------------------------------------------- grouping
    groups = group_by_transcript([v3, v9], ann, "all")
    check("two variants on one transcript are grouped",
          bool(groups) and all(len(v[1]) == 2 for v in groups.values()),
          f"{len(groups)} transcript(s)")
    check("a lone variant is not grouped",
          group_by_transcript([v3], ann, "all") == {})

    # ------------------------------------------------------- the point
    combo = build_combinatorial_proteins([v3, v9], ann, gen,
                                         transcript_mode="all")
    check("a combinatorial protein is produced", bool(combo),
          f"{len(combo)} record(s)")
    if not combo:
        print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
        return 1

    c = combo[0]
    single = {r.sequence for r in (s3 + s9)}
    check("the combined protein is not any single-variant protein",
          c.sequence not in single, f"len={len(c.sequence)}")

    only3 = s3[0].sequence
    only9 = s9[0].sequence
    diff3 = [i for i, (x, y) in enumerate(zip(only3, c.sequence)) if x != y]
    diff9 = [i for i, (x, y) in enumerate(zip(only9, c.sequence)) if x != y]
    check("the combined form differs from each single form at the other "
          "variant's residue", bool(diff3) and bool(diff9),
          f"vs-v3 {diff3[:3]}, vs-v9 {diff9[:3]}")

    check("its class is COMBO", c.variant_class == COMBO_CLASS,
          c.variant_class)
    check("both loci are recorded",
          c.extra["n_variants"] == 2 and len(c.extra["loci"]) == 2,
          str(c.extra.get("loci")))
    check("it is marked unphased, because co-occurrence is not phase",
          "unphased" in c.notes, str(c.notes))
    check("same-class variants are not flagged cross-evidence",
          c.extra["cross_evidence"] is False, str(c.extra["classes"]))

    # ------------------------------------ the case no other tool reaches
    # A somatic SNV and an ADAR edit on one transcript: a proteoform with
    # neither a purely genomic nor a purely transcriptomic basis.
    edit = CombinableVariant("chrT1", p9, r9, a9, "RNA_EDITING", "res_table")
    mixed = build_combinatorial_proteins([v3, edit], ann, gen,
                                         transcript_mode="all")
    check("a DNA variant and an RNA edit combine", bool(mixed))
    if mixed:
        m = mixed[0]
        check("the mix is flagged cross-evidence",
              m.extra["cross_evidence"] is True
              and m.extra["classes"] == ["RNA_EDITING", "SNV"],
              str(m.extra["classes"]))
        check("the evidence types are named in the notes",
              any(n.startswith("cross_evidence_") for n in m.notes),
              str(m.notes))

    # --------------------------------------------------------- refusals
    dup = CombinableVariant("chrT1", p3, r3, a3, "SNV", "vcf")
    over = build_combinatorial_proteins([v3, dup], ann, gen,
                                        transcript_mode="all")
    check("overlapping variants are refused, not mangled", over == [],
          f"{len(over)} record(s)")

    many = []
    for i in range(2, 12):
        gp = gpos_plus(meta, "GPLUS", i, 1)
        b = gen.fetch("chrT1", gp, gp)
        many.append(CombinableVariant("chrT1", gp, b, SWAP[b], "SNV", "vcf"))
    check("a transcript above the variant cap is skipped as an artefact",
          build_combinatorial_proteins(many, ann, gen, transcript_mode="all",
                                       max_variants=3) == [])

    check("an empty variant list yields nothing",
          build_combinatorial_proteins([], ann, gen) == [])

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())

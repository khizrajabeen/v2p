#!/usr/bin/env python3
"""Milestone 5 acceptance test - non-canonical ORFs.

Run:  python tests/test_noncanonical.py
Exit code 0 = all pass. No pytest dependency so it runs anywhere.

The construction is deliberate: a lncRNA carrying a 40-residue ORF in
frame 2 and a 20-residue ORF, so the length floor and the frame bookkeeping
are both checked against a known answer rather than against whatever the
code happens to produce.

The assertion the milestone rests on is the last one - that none of this
appears unless it is asked for. Three-frame translation multiplies database
size, and an inflated search space costs sensitivity.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.build.noncanonical import (                        # noqa: E402
    DEFAULT_MIN_AA, NC_CLASSES, build_noncanonical_proteins, build_utr_orfs,
    classify_biotype, find_orfs,
)
from v2p.seqops import CODON_TABLE, revcomp                  # noqa: E402

PASS: list[str] = []
FAIL: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{' | ' + detail if detail else ''}")


# Reverse codon table, so an ORF can be built from the protein we want.
_AA2CODON: dict[str, str] = {}
for _c, _a in CODON_TABLE.items():
    _AA2CODON.setdefault(_a, _c)


def orf_nt(protein: str) -> str:
    """Nucleotides encoding `protein`, starting ATG, ending TAA."""
    return "ATG" + "".join(_AA2CODON[a] for a in protein[1:]) + "TAA"


@dataclass
class FakeTx:
    """Enough of annotation.Transcript for the builder to work on."""
    tx_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]
    cds: list[tuple[int, int]] = field(default_factory=list)
    tx_biotype: str = "lncRNA"
    gene_name: str = "LINCTEST"
    gene_id: str = "ENSGTEST"


class FakeGenome:
    """Serves one contig from a string, matching Genome.blocks()."""

    def __init__(self, seq: str):
        self.seq = seq

    def blocks(self, chrom, blocks, strand):
        out = "".join(self.seq[s - 1:e] for s, e in sorted(blocks))
        return revcomp(out) if strand == "-" else out


def main() -> int:
    # A 40-residue ORF placed so it begins at offset 2 -> frame 2.
    p40 = "M" + "ACDEFGHIKLMNPQRSTVW"[:19] + "ACDEFGHIKLMNPQRSTVWY"
    assert len(p40) == 40, len(p40)
    p20 = "M" + "ACDEFGHIKLMNPQRSTVW"
    assert len(p20) == 20, len(p20)

    # "GG" pads the front so the 40-mer starts at index 2; a stop separates
    # the two ORFs so they cannot run together.
    seq = "GG" + orf_nt(p40) + "TAA" + orf_nt(p20) + "TAA"
    gen = FakeGenome(seq)
    tx = FakeTx("ENSTNC1", "chrN", "+", [(1, len(seq))])

    # ------------------------------------------------------------- ORFs
    orfs = find_orfs(seq, min_aa=DEFAULT_MIN_AA)
    check("the 40-residue ORF is found", any(o.protein == p40 for o in orfs),
          f"{len(orfs)} orf(s), lengths {[o.length for o in orfs]}")
    hit = next((o for o in orfs if o.protein == p40), None)
    check("its frame is recorded, and it is frame 2",
          hit is not None and hit.frame == 2,
          f"frame={hit.frame if hit else 'n/a'}")
    check("its nucleotide offset points at the ATG",
          hit is not None and seq[hit.nt_start:hit.nt_start + 3] == "ATG",
          seq[hit.nt_start:hit.nt_start + 3] if hit else "n/a")
    check("it is reported as terminating in a stop",
          hit is not None and hit.stop_found)

    check("the 20-residue ORF is rejected at the default threshold",
          all(o.protein != p20 for o in orfs)
          and all(o.length >= DEFAULT_MIN_AA for o in orfs),
          f"lengths {[o.length for o in orfs]}")
    check("lowering the threshold finds the 20-residue ORF",
          any(o.protein == p20 for o in find_orfs(seq, min_aa=20)))
    check("raising the threshold above 40 finds neither",
          find_orfs(seq, min_aa=41) == [])
    lens = [o.length for o in find_orfs(seq, min_aa=5)]
    check("ORFs come back longest first", lens == sorted(lens, reverse=True),
          str(lens))

    # -------------------------------------------------------- start policy
    check("first-ATG policy marks every ORF as ATG-started",
          all(o.started_at_atg for o in orfs))
    no_atg = find_orfs(seq, min_aa=DEFAULT_MIN_AA, require_atg=False)
    check("any-start policy records the frame and no ATG claim",
          no_atg and all(not o.started_at_atg for o in no_atg)
          and all(o.frame in (0, 1, 2) for o in no_atg),
          f"{len(no_atg)} orf(s)")

    # ---------------------------------------------------------- biotypes
    check("lncRNA maps to NC_LNCRNA",
          classify_biotype("lncRNA") == "NC_LNCRNA")
    check("a pseudogene maps to NC_PSEUDOGENE",
          classify_biotype("processed_pseudogene") == "NC_PSEUDOGENE"
          and classify_biotype("transcribed_unitary_pseudogene")
          == "NC_PSEUDOGENE")
    check("an unrecognised non-coding biotype falls back to NC_ALTORF",
          classify_biotype("some_new_gencode_biotype") == "NC_ALTORF")
    check("protein_coding is skipped",
          classify_biotype("protein_coding") is None)
    check("structural RNAs are skipped",
          classify_biotype("miRNA") is None
          and classify_biotype("snoRNA") is None
          and classify_biotype("rRNA") is None)
    check("an empty biotype is skipped", classify_biotype("") is None)

    # ------------------------------------------------------------ records
    recs = build_noncanonical_proteins([tx], gen)
    check("one record is built for the transcript", len(recs) == 1,
          f"{len(recs)} record(s)")
    r = recs[0]
    check("the record carries the 40-residue ORF", r.sequence == p40,
          f"len={len(r.sequence)}")
    check("its class is NC_LNCRNA and in the documented vocabulary",
          r.variant_class == "NC_LNCRNA" and r.variant_class in NC_CLASSES,
          r.variant_class)
    check("the frame survives into the record",
          r.extra.get("frame") == 2 and "frame2" in r.notes,
          f"{r.extra.get('frame')}/{r.notes}")
    check("the consequence names it a non-canonical ORF",
          r.consequence == "noncanonical_orf", r.consequence)
    check("the transcript id is carried", r.transcript == "ENSTNC1")

    # A coding transcript must never be three-frame translated: its real
    # protein is already emitted by the normal path.
    coding = FakeTx("ENSTC1", "chrN", "+", [(1, len(seq))],
                    cds=[(3, 123)], tx_biotype="protein_coding")
    check("a coding transcript is not three-frame translated",
          build_noncanonical_proteins([coding], gen) == [])
    check("max_per_transcript caps the output",
          len(build_noncanonical_proteins([tx], gen, min_aa=5,
                                          max_per_transcript=2)) == 2)

    # ------------------------------------------------ minus strand sanity
    rc = revcomp(seq)
    txm = FakeTx("ENSTNC2", "chrN", "-", [(1, len(rc))])
    check("a minus-strand transcript finds the same ORF",
          any(x.sequence == p40
              for x in build_noncanonical_proteins([txm], FakeGenome(rc))))

    # ----------------------------------------------------------- NC_UTR
    # uORFs in the 5-prime UTR of a CODING transcript. The class existed in
    # the vocabulary before anything emitted it; this is that gap closed.
    class UtrTx(FakeTx):
        def cds_offset_in_tx(self):
            return self._cds_off

    utr_nt = orf_nt(p40)                       # the uORF, 123 nt
    cds_nt = orf_nt("M" + "ACDEFGHIK" * 3)     # the real CDS
    full = utr_nt + cds_nt + orf_nt(p40)       # 5-prime UTR, CDS, 3-prime UTR
    utx = UtrTx("ENSTU1", "chrN", "+", [(1, len(full))],
                cds=[(len(utr_nt) + 1, len(utr_nt) + len(cds_nt))],
                tx_biotype="protein_coding")
    utx._cds_off = len(utr_nt)
    urecs = build_utr_orfs([utx], FakeGenome(full))
    check("NC_UTR records are produced for a coding transcript",
          len(urecs) >= 1, f"{len(urecs)} record(s)")
    check("every NC_UTR record is classed NC_UTR",
          all(r.variant_class == "NC_UTR" for r in urecs))
    check("the 5-prime uORF is found and its region recorded",
          any(r.sequence == p40 and r.extra.get("utr") == "5_prime"
              for r in urecs),
          str([(len(r.sequence), r.extra.get("utr")) for r in urecs]))
    check("the UTR region is named in the notes",
          all(any(n.endswith("_utr_orf") for n in r.notes) for r in urecs))
    check("the main CDS is not re-emitted as a UTR ORF",
          all("ACDEFGHIKACDEFGHIK" not in r.sequence for r in urecs))
    check("a non-coding transcript yields no UTR ORFs",
          build_utr_orfs([tx], gen) == [])

    # ------------------------------------------------------ off by default
    from v2p.fasta import VARIANT_TYPE_VOCAB                 # noqa: E402
    check("NC classes are absent from the default variant vocabulary",
          not any(c in VARIANT_TYPE_VOCAB for c in NC_CLASSES),
          str([c for c in NC_CLASSES if c in VARIANT_TYPE_VOCAB]))

    stage = (ROOT / "src" / "v2p" / "stages"
             / "02_build_protein_fasta.py").read_text(encoding="utf-8")
    check("the stage exposes --include-noncanonical",
          "--include-noncanonical" in stage)
    if "--include-noncanonical" in stage:
        after = stage.split("--include-noncanonical")[1][:220]
        check("the flag is store_true, so it defaults off",
              "store_true" in after and "default=True" not in after,
              after.strip()[:70])

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())

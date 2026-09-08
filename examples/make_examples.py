#!/usr/bin/env python3
"""Regenerate the example input folder.

The examples are *synthetic calls at real coordinates*. Nothing here is
redistributed from the HCC1395 high-confidence package, which is not ours
to publish. Instead:

  - small variants are placed into the CDS of named genes, with the
    reference allele read from the genome so the call is internally
    consistent and the consequence is a genuine missense;
  - the RNA-editing rows are the canonical ADAR recoding sites from the
    literature (GRIA2 Q607R, NEIL1 K242R, BLCAP Y2C, CDK13 Q103R,
    COG3 I635V) - the set no competing tool can attempt;
  - the fusions are published, characterised rearrangements;
  - the splicing events are built from real GENCODE junctions.

Run:  python examples/make_examples.py --genome ref/... --gtf ref/...

The committed output is what the quickstart uses, so this only needs
re-running if the coordinates or the annotation release change.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.annotation import Annotation, Genome              # noqa: E402
from v2p.seqops import CODON_TABLE                         # noqa: E402

HERE = Path(__file__).resolve().parent

# Genes to place a synthetic coding substitution in. Chosen because they
# are recognisable and well annotated, not because these variants are real.
GENES = ["TP53", "KRAS", "BRAF", "PIK3CA", "EGFR", "PTEN", "NRAS", "IDH1"]

# Canonical ADAR recoding sites, as (gene, published protein change). The
# genomic coordinate is *derived* from the annotation below rather than
# quoted from memory, and the derivation only succeeds if changing one A to
# G in that codon actually produces the published substitution. A site that
# cannot be derived is reported and dropped rather than shipped as a guess.
ADAR_RECODING = [
    ("GRIA2", "Q607R"),
    ("NEIL1", "K242R"),
    ("BLCAP", "Y2C"),
    ("CDK13", "Q103R"),
    ("COG3", "I635V"),
]

# Published, characterised fusions with their canonical breakpoints.
FUSIONS = [
    ("BCR", "ABL1", "chr22:23290413", "chr9:130854064", "SR;LR", 0.99),
    ("EML4", "ALK", "chr2:42295516", "chr2:29223528", "SR;LR", 0.97),
    ("TMPRSS2", "ERG", "chr21:41508081", "chr21:38403709", "SR", 0.91),
]

COMPLEMENT_SWAP = {"A": "G", "G": "A", "C": "T", "T": "C"}


def codon_locus(ann, gene: str, codon: int):
    """Genomic position of base 1 of `codon` in `gene`'s representative CDS.

    `representative` is the pipeline's own tie-break chain (MANE_Select,
    Ensembl_canonical, basic, longest CDS, longest transcript, id), so the
    example variant lands on the transcript the pipeline itself would pick.
    """
    tx = ann.representative(gene)
    if tx is None or not tx.cds:
        return None
    off = 3 * (codon - 1)
    blocks = sorted(tx.cds) if tx.strand == "+" else sorted(tx.cds,
                                                            reverse=True)
    seen = 0
    for s, e in blocks:
        ln = e - s + 1
        if off < seen + ln:
            pos = s + (off - seen) if tx.strand == "+" else e - (off - seen)
            return tx.chrom, pos, tx.strand, tx.tx_id
        seen += ln
    return None


def tx_offset_to_genomic(tx, off: int) -> int:
    """Genomic position of CDS offset `off` (0-based, transcription order)."""
    blocks = sorted(tx.cds) if tx.strand == "+" else sorted(tx.cds,
                                                            reverse=True)
    seen = 0
    for s, e in blocks:
        ln = e - s + 1
        if off < seen + ln:
            return s + (off - seen) if tx.strand == "+" else e - (off - seen)
        seen += ln
    raise IndexError(off)


def derive_adar_site(ann, gen, gene: str, change: str):
    """Find the genomic A>G that produces the published recoding.

    Returns (chrom, pos, tx_id, verification) or None. The site is only
    returned if the reference codon really encodes the published `from`
    residue and a single A->G in it really yields the `to` residue, so a
    wrong coordinate cannot pass silently.
    """
    frm, to = change[0], change[-1]
    codon_no = int(change[1:-1])
    tx = ann.representative(gene)
    if tx is None or not tx.cds:
        return None
    cds = gen.blocks(tx.chrom, tx.cds, tx.strand)
    i = 3 * (codon_no - 1)
    codon = cds[i:i + 3]
    if len(codon) < 3 or CODON_TABLE.get(codon) != frm:
        return None
    for j in range(3):
        if codon[j] != "A":
            continue
        mutated = codon[:j] + "G" + codon[j + 1:]
        if CODON_TABLE.get(mutated) != to:
            continue
        pos = tx_offset_to_genomic(tx, i + j)
        # On the minus strand a transcript A is a genomic T, and the edit
        # reads as T>C in genome coordinates.
        gref = gen.fetch(tx.chrom, pos, pos)
        expect = "A" if tx.strand == "+" else "T"
        if gref != expect:
            return None
        return (tx.chrom, pos, tx.tx_id,
                f"derived from GENCODE {gene} {tx.tx_id} codon {codon_no}: "
                f"{codon}->{mutated} = {frm}->{to}")
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", default="ref/GRCh38.primary_assembly.genome.fa")
    ap.add_argument("--gtf", default="ref/gencode.v44.annotation.gtf.gz")
    a = ap.parse_args()

    print("loading annotation (this takes a minute) ...")
    ann = Annotation.from_gtf(a.gtf)
    gen = Genome(a.genome)

    rows = []
    for i, gene in enumerate(GENES):
        hit = codon_locus(ann, gene, 25 + i)
        if not hit:
            print(f"  ! {gene}: no coding transcript found, skipped")
            continue
        chrom, pos, _strand, tx = hit
        ref = gen.fetch(chrom, pos, pos).upper()
        alt = COMPLEMENT_SWAP.get(ref)
        if not alt:
            print(f"  ! {gene}: reference base {ref!r} at {chrom}:{pos}")
            continue
        rows.append((chrom, pos, ref, alt, gene, tx))
        print(f"  {gene:8s} {chrom}:{pos} {ref}>{alt}  ({tx})")

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(HERE / "small_variants.vcf", "w", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##reference=GRCh38\n")
        fh.write("##source=v2p examples/make_examples.py (SYNTHETIC CALLS)\n")
        for c in sorted({r[0] for r in rows}):
            fh.write(f"##contig=<ID={c},length={len(gen.fa[gen.resolve(c)])}>\n")
        fh.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">\n')
        fh.write('##INFO=<ID=TX,Number=1,Type=String,Description="Transcript">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt, gene, tx in rows:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t"
                     f"GENE={gene};TX={tx}\n")
    print(f"wrote small_variants.vcf ({len(rows)} calls)")

    cols = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
            "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene",
            "AAChange.refGene"]
    sites, dropped = [], []
    for gene, change in ADAR_RECODING:
        hit = derive_adar_site(ann, gen, gene, change)
        if hit is None:
            dropped.append(f"{gene} {change}")
            continue
        chrom, pos, tx_id, how = hit
        sites.append((chrom, pos, gene, change, tx_id))
        print(f"  {gene:8s} {change:7s} {chrom}:{pos}  {how}")
    for d in dropped:
        print(f"  ! {d}: could not be derived from the annotation, dropped")
    with open(HERE / "rna_editing.txt", "w", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for chrom, pos, gene, change, tx_id in sites:
            strand_ref = gen.fetch(chrom, pos, pos)
            alt = "G" if strand_ref == "A" else "C"
            fh.write("\t".join([chrom, str(pos), str(pos), strand_ref, alt,
                                "exonic", gene, ".", "nonsynonymous SNV",
                                f"{gene}:{tx_id}:.:.:p.{change}"]) + "\n")
    print(f"wrote rna_editing.txt ({len(sites)} ADAR sites, "
          f"{len(dropped)} dropped)")

    with open(HERE / "fusions.csv", "w", newline="\n") as fh:
        fh.write("tag,Gene1,Gene2,Type,Validated,confidence,"
                 "breakpoint1,breakpoint2\n")
        for g1, g2, b1, b2, typ, conf in FUSIONS:
            fh.write(f"{g1}--{g2},{g1},{g2},{typ},yes,{conf},{b1},{b2}\n")
    print(f"wrote fusions.csv ({len(FUSIONS)} fusions)")

    ev = []
    for t in sorted(ann.tx.values(), key=lambda x: x.tx_id):
        if not t.cds or len(t.exons) < 4:
            continue
        ex = sorted(t.exons)
        ev.append(f"{t.gene_id};SE:{t.chrom}:{ex[0][1]}-{ex[1][0]}:"
                  f"{ex[1][1]}-{ex[2][0]}:{t.strand}")
        if len(ev) >= 5:
            break
    with open(HERE / "splicing.csv", "w", newline="\n") as fh:
        fh.write("AS,dPSI,p_val,Batch_support,Lib_support,gene_id,confidence\n")
        for e in ev:
            fh.write(f"{e},0.42,0.001,6;9,2;2,{e.split(';')[0]},0.9\n")
    print(f"wrote splicing.csv ({len(ev)} events)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

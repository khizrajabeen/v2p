#!/usr/bin/env python3
"""Build a miniature genome + GTF with analytically known proteins.

Two genes are laid down, one per strand, each with two exons, plus a
third single-exon gene used as a fusion 3' partner. Because the
transcript sequence is constructed first and then *placed* into the
contigs, the expected protein for every test case is known by
construction rather than by trusting the code under test.
"""

from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIX = HERE / "fixtures"

# one unambiguous codon per amino acid
CODON = {
    "U": "TGA",
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT",
    "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT",
    "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
    "W": "TGG", "Y": "TAT", "*": "TAA",
}

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def cds_for(protein: str) -> str:
    return "".join(CODON[a] for a in protein) + CODON["*"]


def filler(n: int, rng: random.Random) -> str:
    """Intergenic/intronic filler with no ATG, so stray ORFs can't appear."""
    out = []
    while len(out) < n:
        b = rng.choice("ACGT")
        if len(out) >= 2 and out[-2] == "A" and out[-1] == "T" and b == "G":
            continue
        out.append(b)
    return "".join(out)


def build() -> dict:
    rng = random.Random(20260903)
    FIX.mkdir(parents=True, exist_ok=True)

    # ---- gene definitions ------------------------------------------------
    prot_plus = "MASTKLYVWQDEGHRNFIPC"
    prot_minus = "MGKTVLYSEDAWQRNFHIP"
    prot_part3 = "MQVLKDEFGHTYWNRSAIP"
    # GSEC: a selenoprotein - residue 6 is encoded by TGA, which GENCODE
    # annotates as Selenocysteine. GMITO: a chrM gene using table 2.
    prot_sec = "MAKDLUGHTYWNRSAIPQV"
    prot_mito = "MWKTIVLYSEDAWQRNFHI"

    specs = [
        # name, chrom, strand, protein, utr5, utr3, exon1 start, intron len
        ("GPLUS", "chrT1", "+", prot_plus, 21, 30, 101, 200),
        ("GMINUS", "chrT2", "-", prot_minus, 21, 30, 101, 200),
        ("GPART3", "chrT3", "+", prot_part3, 21, 30, 101, 200),
        ("GSEC", "chrT4", "+", prot_sec, 21, 30, 101, 200),
        ("GMITO", "chrM", "+", prot_mito, 21, 30, 101, 200),
        # GDUP is deliberately declared twice, on chrT5 and chrT6, to
        # reproduce an ambiguous gene symbol.
        ("GDUP", "chrT5", "+", prot_plus, 21, 30, 101, 200),
        ("GDUP", "chrT6", "-", prot_minus, 21, 30, 101, 200),
    ]

    contigs: dict[str, list[str]] = {}
    gtf_lines: list[str] = []
    meta: dict = {}

    for name, chrom, strand, prot, u5, u3, ex1_start, intron in specs:
        cds = cds_for(prot)
        tx = filler(u5, rng) + cds + filler(u3, rng)
        # split so the junction falls inside the CDS, not on a codon boundary
        split = u5 + 3 * (len(prot) // 2) + 1
        b1, b2 = tx[:split], tx[split:]

        contig = list(filler(3000, rng))

        if strand == "+":
            e1s, e1e = ex1_start, ex1_start + len(b1) - 1
            e2s, e2e = e1e + intron + 1, e1e + intron + len(b2)
            for i, ch in enumerate(b1):
                contig[e1s - 1 + i] = ch
            for i, ch in enumerate(b2):
                contig[e2s - 1 + i] = ch
            cds_gstart = e1s + u5
        else:
            # transcript runs right-to-left: block1 is the rightmost exon
            e2s, e2e = ex1_start, ex1_start + len(b2) - 1          # 3' exon
            e1s, e1e = e2e + intron + 1, e2e + intron + len(b1)    # 5' exon
            for i, ch in enumerate(revcomp(b1)):
                contig[e1s - 1 + i] = ch
            for i, ch in enumerate(revcomp(b2)):
                contig[e2s - 1 + i] = ch
            cds_gstart = e1e - u5

        contigs[chrom] = contig
        exons = sorted([(e1s, e1e), (e2s, e2e)])
        # CDS blocks
        cds_blocks = []
        if strand == "+":
            cds_end_tx = u5 + len(cds)
            remaining = len(cds)
            pos = cds_gstart
            for s, e in exons:
                if pos > e:
                    continue
                take = min(e - pos + 1, remaining)
                cds_blocks.append((pos, pos + take - 1))
                remaining -= take
                if remaining == 0:
                    break
                pos = exons[1][0]
            _ = cds_end_tx
        else:
            remaining = len(cds)
            pos = cds_gstart
            for s, e in sorted(exons, reverse=True):
                if pos < s:
                    continue
                take = min(pos - s + 1, remaining)
                cds_blocks.append((pos - take + 1, pos))
                remaining -= take
                if remaining == 0:
                    break
                pos = exons[0][1]

        gid, tid = f"{name}_G_{chrom}", f"{name}_T_{chrom}"
        attrs = (f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{name}"; '
                 f'gene_type "protein_coding"; transcript_type "protein_coding"; '
                 f'tag "basic"; tag "Ensembl_canonical";')
        for s, e in exons:
            gtf_lines.append(f"{chrom}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t{attrs}")
        for s, e in sorted(cds_blocks):
            gtf_lines.append(f"{chrom}\ttest\tCDS\t{s}\t{e}\t0\t{strand}\t0\t{attrs}")
        if "U" in prot:
            i = prot.index("U")
            g = cds_gstart + 3 * i          # + strand fixtures only
            gtf_lines.append(f"{chrom}\ttest\tSelenocysteine\t{g}\t{g+2}"
                             f"\t.\t{strand}\t.\t{attrs}")

        meta[f"{name}@{chrom}"] = {
            "chrom": chrom, "strand": strand, "protein": prot,
            "exons": exons, "cds_blocks": sorted(cds_blocks),
            "utr5": u5, "cds_gstart": cds_gstart, "tx": tx, "cds": cds,
        }

    fa = FIX / "mini.fa"
    with open(fa, "w") as fh:
        for chrom, seq in contigs.items():
            fh.write(f">{chrom}\n")
            s = "".join(seq)
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    (FIX / "mini.gtf").write_text("\n".join(gtf_lines) + "\n")
    return meta


if __name__ == "__main__":
    m = build()
    for k, v in m.items():
        print(k, v["chrom"], v["strand"], "protein:", v["protein"])
        print("   exons", v["exons"], "cds", v["cds_blocks"])

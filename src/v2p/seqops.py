"""Nucleotide -> amino-acid machinery.

Deliberately dependency-light: the standard genetic code is hard-coded so
the translation step is auditable and does not silently change if a
library updates. Biopython is used only where it adds value elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass

# NCBI translation table 1 (Standard). Selenocysteine / pyrrolysine
# recoding is NOT applied: stop codons terminate translation. Any
# selenoprotein in the input will be truncated at the UGA; such
# transcripts are flagged rather than silently mistranslated.
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# NCBI translation table 2 (Vertebrate Mitochondrial). Differences from
# table 1: TGA is tryptophan rather than a stop, ATA is methionine rather
# than isoleucine, and AGA/AGG are stops rather than arginine. chrM
# transcripts translated with table 1 read straight through their real
# stops and past their real tryptophans.
CODON_TABLE_MITO: dict[str, str] = dict(CODON_TABLE)
CODON_TABLE_MITO.update({"TGA": "W", "ATA": "M", "AGA": "*", "AGG": "*"})

MITO_CONTIGS = {"chrM", "chrMT", "MT", "M"}


def table_for_contig(chrom: str) -> dict[str, str]:
    """The genetic code appropriate to a contig."""
    return CODON_TABLE_MITO if chrom in MITO_CONTIGS else CODON_TABLE


_COMPLEMENT = str.maketrans("ACGTNacgtnRYKMBVDHrykmbvdh",
                            "TGCANtgcanYRMKVBHDyrmkvbhd")

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
}


def revcomp(seq: str) -> str:
    """Reverse complement, IUPAC-aware."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class Translation:
    """Result of translating an ORF."""

    protein: str          # peptide up to (not including) the first stop
    stop_found: bool      # True if a stop codon terminated translation
    trailing_partial: int  # nt left over that did not form a full codon
    ambiguous_codons: int  # codons containing non-ACGT characters -> 'X'
    selenocysteines: int = 0   # UGA codons recoded to U


def translate_orf(nt: str, table: dict[str, str] | None = None,
                  to_stop: bool = True,
                  sec_codons: set[int] | None = None) -> Translation:
    """Translate `nt` in frame 0.

    Codons containing any non-ACGT base become 'X' rather than raising:
    reference assemblies contain N runs and refusing to translate them
    would silently drop real transcripts.

    `sec_codons` holds 0-based codon indices that GENCODE annotates as
    selenocysteine. A UGA at one of those positions is recoded to U
    instead of terminating translation. Without this, all 25
    selenoprotein genes are truncated at their first Sec codon, which is
    not a subtle error: SELENOP alone loses most of its length.
    """
    tab = table or CODON_TABLE
    sec = sec_codons or set()
    nt = nt.upper().replace("U", "T")
    aas: list[str] = []
    ambiguous = 0
    n_sec = 0
    stop_found = False
    n_full = len(nt) // 3
    for i in range(n_full):
        codon = nt[3 * i: 3 * i + 3]
        aa = tab.get(codon)
        if aa is None:
            aa = "X"
            ambiguous += 1
        if aa == "*":
            if i in sec and codon == "TGA":
                aa = "U"
                n_sec += 1
            else:
                stop_found = True
                if to_stop:
                    break
        aas.append(aa)
    return Translation(
        protein="".join(aas),
        stop_found=stop_found,
        trailing_partial=len(nt) - 3 * n_full,
        ambiguous_codons=ambiguous,
        selenocysteines=n_sec,
    )


def first_atg_offset(nt: str, min_aa: int = 8) -> int | None:
    """Index of the first ATG that opens an ORF of >= `min_aa` residues.

    Used for isoforms whose annotated start codon is destroyed by the
    variant (e.g. alternative-first-exon events).
    """
    seq = nt.upper()
    for i in range(len(seq) - 2):
        if seq[i:i + 3] != "ATG":
            continue
        if len(translate_orf(seq[i:]).protein) >= min_aa:
            return i
    return None


def peptide_window(protein: str, pos0: int, flank: int) -> tuple[str, int]:
    """`flank` residues either side of 0-based `pos0`.

    Returns (subsequence, index-of-pos0-within-subsequence). Used to emit
    short variant-centred peptides instead of whole proteins when the
    downstream search engine wants a compact database.
    """
    start = max(0, pos0 - flank)
    end = min(len(protein), pos0 + flank + 1)
    return protein[start:end], pos0 - start


def first_difference(ref: str, alt: str) -> int:
    """0-based index of the first differing residue, else min length."""
    n = min(len(ref), len(alt))
    for i in range(n):
        if ref[i] != alt[i]:
            return i
    return n


def hgvs_p_to_change(hgvs: str) -> tuple[str, int, str] | None:
    """Parse a simple ANNOVAR/HGVS protein change: 'p.K509E' -> ('K',509,'E').

    Returns None for anything that is not a simple substitution (frameshift,
    duplication, delins), which the caller must handle by full re-translation.
    """
    s = hgvs.strip()
    if s.startswith("p."):
        s = s[2:]
    if not s:
        return None
    # three-letter form, e.g. Lys509Glu
    for three, one in AA3_TO_1.items():
        if s.startswith(three):
            rest = s[len(three):]
            digits = ""
            for ch in rest:
                if ch.isdigit():
                    digits += ch
                else:
                    break
            tail = rest[len(digits):]
            if digits and tail in AA3_TO_1:
                return one, int(digits), AA3_TO_1[tail]
    # one-letter form, e.g. K509E
    if len(s) >= 3 and s[0].isalpha() and s[-1] in "ACDEFGHIKLMNPQRSTVWY*X":
        digits = s[1:-1]
        if digits.isdigit():
            return s[0], int(digits), s[-1]
    return None

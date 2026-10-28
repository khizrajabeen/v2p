"""Peptide-level operations for search-database preparation.

Three things the proteogenomics literature treats as mandatory and that
a protein-level pipeline does not give you for free:

1. **Peptide-level novelty.** A variant protein entry is only useful if
   it contributes at least one tryptic peptide absent from the reference
   proteome. An entry whose every peptide already exists adds search
   space and no information. Reporting this per entry lets the recipient
   filter on evidence rather than on faith.

2. **Decoys.** FDR control needs a decoy database built by the same
   process as the target. Reversing while holding the C-terminal residue
   fixed ("pseudo-reverse") preserves the tryptic peptide mass
   distribution, which reversed-whole-protein decoys do not when the
   database is peptide-derived.

3. **Class-specific FDR support.** Variant peptides need their own FDR
   estimation, separate from canonical peptides. That requires the class
   label to survive into the search output, which is what the `VT=` field
   in every header is for.
"""

from __future__ import annotations

import re
from collections import Counter

# Trypsin: cleave after K or R, not before P.
_TRYPSIN = re.compile(r"(?<=[KR])(?!P)")


def digest(seq: str, missed_cleavages: int = 2,
           min_len: int = 6, max_len: int = 50) -> set[str]:
    """In-silico tryptic peptides, allowing missed cleavages."""
    seq = seq.replace("*", "")
    if not seq:
        return set()
    frags = [f for f in _TRYPSIN.split(seq) if f]
    out: set[str] = set()
    for i in range(len(frags)):
        pep = ""
        for j in range(i, min(i + missed_cleavages + 1, len(frags))):
            pep += frags[j]
            if min_len <= len(pep) <= max_len:
                out.add(pep)
            elif len(pep) > max_len:
                break
    return out


def reference_peptide_space(sequences, missed_cleavages: int = 2,
                            min_len: int = 6, max_len: int = 50,
                            equate_il: bool = True) -> set[str]:
    """All tryptic peptides of the reference proteome.

    `equate_il` folds isoleucine and leucine together. They are
    isobaric, so a mass spectrometer cannot distinguish them and a
    peptide differing only by I/L is not detectably novel. Treating them
    as distinct overstates novelty.
    """
    space: set[str] = set()
    for s in sequences:
        for p in digest(s, missed_cleavages, min_len, max_len):
            space.add(p.replace("I", "L") if equate_il else p)
    return space


def novel_peptides(seq: str, reference_space: set[str],
                   missed_cleavages: int = 2, min_len: int = 6,
                   max_len: int = 50, equate_il: bool = True) -> set[str]:
    """Tryptic peptides of `seq` that are absent from the reference."""
    out = set()
    for p in digest(seq, missed_cleavages, min_len, max_len):
        key = p.replace("I", "L") if equate_il else p
        if key not in reference_space:
            out.add(p)
    return out


def variant_spanning_peptides(seq: str, pos_aa: int | None,
                              novel: set[str]) -> set[str]:
    """Novel peptides that actually cover the variant residue.

    A frameshift protein shares its whole N-terminal region with the
    reference; only the peptides overlapping or downstream of the
    variant carry evidence for it. These are the peptides worth reporting
    as identifying the variant.
    """
    if pos_aa is None:
        return novel
    out = set()
    for p in novel:
        i = seq.find(p)
        while i != -1:
            if i < pos_aa <= i + len(p):
                out.add(p)
                break
            i = seq.find(p, i + 1)
    return out


# --------------------------------------------------------------------------
# decoys
# --------------------------------------------------------------------------

def decoy_reverse(seq: str) -> str:
    """Plain reversal of the whole sequence."""
    return seq[::-1]


def decoy_pseudo_reverse(seq: str) -> str:
    """Reverse each tryptic peptide while keeping its C-terminal K/R.

    This preserves the number and mass distribution of tryptic peptides,
    so the decoy population models false positives from the same peptide
    space as the target. A whole-protein reversal changes the cleavage
    pattern and biases FDR estimation for variant-peptide subsets.
    """
    frags = [f for f in _TRYPSIN.split(seq) if f]
    out = []
    for f in frags:
        if f and f[-1] in "KR":
            out.append(f[:-1][::-1] + f[-1])
        else:
            out.append(f[::-1])
    return "".join(out)


def decoy_shuffle(seq: str, seed: int = 0, max_tries: int = 10) -> str:
    """Deterministic within-peptide shuffle, C-terminal residue fixed.

    Every fragment is shuffled, including the protein's last one, which
    ends at no cleavage site: skipping it copied a genuine target peptide
    into the decoy once per protein, and target/decoy overlap biases the
    FDR these decoys exist to estimate.

    Proline holds its position too. Trypsin does not cleave K-P or R-P,
    so a free proline can move next to a K and destroy a cleavage site,
    leaving the decoy with a different peptide count from the target.
    OpenMS pins K, R and P for the same reason.

    A shuffle can land back on the input, so it is retried. A fragment
    with fewer than two free positions is returned unchanged.
    """
    import random

    rng = random.Random(seed)
    out = []
    for f in (f for f in _TRYPSIN.split(seq) if f):
        body, tail = f[:-1], f[-1]
        free = [i for i, aa in enumerate(body) if aa != "P"]
        if len(free) < 2:
            out.append(f)
            continue
        shuffled = f
        for _ in range(max_tries):
            order = free[:]
            rng.shuffle(order)
            chars = list(body)
            for dst, src in zip(free, order):
                chars[dst] = body[src]
            shuffled = "".join(chars) + tail
            if shuffled != f:
                break
        out.append(shuffled)
    return "".join(out)


DECOY_METHODS = {
    "reverse": decoy_reverse,
    "pseudo_reverse": decoy_pseudo_reverse,
    "shuffle": decoy_shuffle,
}


def digest_stats(sequences) -> dict:
    """Peptide-space summary, for reporting database size honestly."""
    peps: Counter = Counter()
    for s in sequences:
        for p in digest(s):
            peps[len(p)] += 1
    total = sum(peps.values())
    return {"n_peptides": total,
            "mean_length": round(sum(k * v for k, v in peps.items())
                                 / total, 2) if total else 0}

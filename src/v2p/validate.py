"""UniProt reference handling and independent validation of translations.

The UniProt human canonical proteome is used three ways:

  1. **Format template.** The delivered FASTA must match its header
     grammar and its single-line sequence layout.
  2. **Independent check on the annotation.** For every variant with a
     reported protein change, the *reference* residue at that position is
     compared against UniProt. A mismatch means the annotation's
     transcript and UniProt's canonical isoform disagree, or the call is
     wrong. Either way it must be seen, not assumed away.
  3. **Ground truth for the translation engine.** Reference proteins
     produced by translating GENCODE transcripts are compared residue by
     residue against UniProt. If our GRIA2 reference protein does not
     match UniProt's GRIA2, the variant GRIA2 protein is not trustworthy
     either. This is the single most informative QC gate in the pipeline
     and it needs no orthogonal tool to run.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

_HDR = re.compile(
    r"^>(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<desc>.*?)"
    r"(?=\s+(?:OS|OX|GN|PE|SV)=|$)"
)
_KV = re.compile(r"\b(OS|OX|GN|PE|SV)=(.*?)(?=\s+(?:OS|OX|GN|PE|SV)=|$)")


@dataclass
class UniProtEntry:
    db: str
    accession: str
    entry_name: str
    description: str
    os: str = ""
    ox: str = ""
    gn: str = ""
    pe: str = ""
    sv: str = ""
    sequence: str = ""

    def header(self) -> str:
        h = f">{self.db}|{self.accession}|{self.entry_name} {self.description}"
        for k, v in (("OS", self.os), ("OX", self.ox), ("GN", self.gn),
                     ("PE", self.pe), ("SV", self.sv)):
            if v:
                h += f" {k}={v}"
        return h


def parse_uniprot_header(line: str) -> UniProtEntry | None:
    m = _HDR.match(line.rstrip("\n"))
    if not m:
        return None
    e = UniProtEntry(db=m.group("db"), accession=m.group("acc"),
                     entry_name=m.group("entry"),
                     description=m.group("desc").strip())
    for k, v in _KV.findall(line):
        setattr(e, k.lower(), v.strip())
    return e


class UniProtDB:
    """Indexed UniProt proteome."""

    def __init__(self) -> None:
        self.entries: dict[str, UniProtEntry] = {}
        self.by_gene: dict[str, list[str]] = {}
        self.malformed_headers: list[str] = []

    @classmethod
    def load(cls, path: str | Path, logger=None) -> "UniProtDB":
        db = cls()
        p = Path(path)
        opener = gzip.open if p.suffix == ".gz" else open
        cur: UniProtEntry | None = None
        chunks: list[str] = []
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur is not None:
                        cur.sequence = "".join(chunks)
                        db._add(cur)
                    chunks = []
                    cur = parse_uniprot_header(line)
                    if cur is None:
                        db.malformed_headers.append(line.rstrip("\n")[:120])
                elif cur is not None:
                    chunks.append(line.strip())
        if cur is not None:
            cur.sequence = "".join(chunks)
            db._add(cur)
        if logger:
            logger.info("UniProt %s: %d entries, %d genes, %d malformed headers",
                        p.name, len(db.entries), len(db.by_gene),
                        len(db.malformed_headers))
        return db

    def _add(self, e: UniProtEntry) -> None:
        self.entries[e.accession] = e
        if e.gn:
            for g in e.gn.split():
                self.by_gene.setdefault(g, []).append(e.accession)

    def for_gene(self, gene: str) -> list[UniProtEntry]:
        return [self.entries[a] for a in self.by_gene.get(gene, [])]

    def longest_for_gene(self, gene: str) -> UniProtEntry | None:
        c = self.for_gene(gene)
        if not c:
            return None
        return max(c, key=lambda e: (len(e.sequence), e.accession))


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------

@dataclass
class ResidueCheck:
    gene: str
    position: int
    expected_aa: str
    observed_aa: str
    accession: str
    status: str            # match | mismatch | out_of_range | no_uniprot_entry
    uniprot_length: int = 0
    extra: dict = field(default_factory=dict)


def check_reference_residue(db: UniProtDB, gene: str, pos: int,
                            expected_aa: str) -> ResidueCheck:
    """Is `expected_aa` really at 1-based `pos` of this gene's protein?

    All UniProt entries for the gene are tried; a match against any of them
    counts, because the variant caller's transcript may correspond to a
    different canonical entry than the longest one.
    """
    cands = db.for_gene(gene)
    if not cands:
        return ResidueCheck(gene, pos, expected_aa, "", "", "no_uniprot_entry")
    best = ResidueCheck(gene, pos, expected_aa, "", "", "out_of_range")
    for e in cands:
        if pos < 1 or pos > len(e.sequence):
            if best.status == "out_of_range":
                best = ResidueCheck(gene, pos, expected_aa, "", e.accession,
                                    "out_of_range", len(e.sequence))
            continue
        obs = e.sequence[pos - 1]
        rc = ResidueCheck(gene, pos, expected_aa, obs, e.accession,
                          "match" if obs == expected_aa else "mismatch",
                          len(e.sequence))
        if rc.status == "match":
            return rc
        best = rc
    return best


def _anchor_offset_k(ours: str, theirs: str, k: int) -> int | None:
    """Modal offset of `ours` within `theirs` using unique k-mer anchors."""
    if len(ours) < k or len(theirs) < k:
        return None
    n_probe = 12
    span = max(1, (len(ours) - k) // n_probe)
    votes: dict[int, int] = {}
    for start in range(0, len(ours) - k + 1, span):
        kmer = ours[start:start + k]
        j = theirs.find(kmer)
        if j < 0 or theirs.find(kmer, j + 1) != -1:
            continue                      # absent, or not unique
        votes[j - start] = votes.get(j - start, 0) + 1
    if not votes:
        return None
    off, n = max(votes.items(), key=lambda kv: (kv[1], -abs(kv[0])))
    return off if (n >= 2 or len(votes) == 1) else None


def _anchor_offset(ours: str, theirs: str) -> int | None:
    """Offset of `ours` within `theirs`.

    A positional comparison is worthless when two isoforms differ only by
    an N-terminal extension: every residue after the offset counts as a
    mismatch and identity collapses toward zero even though the sequences
    are the same protein. Anchoring on interior k-mers recovers the real
    relationship in O(n). Shorter anchors are tried in turn because a
    densely substituted sequence may leave no clean 30-mer.
    """
    for k in (30, 18, 12):
        off = _anchor_offset_k(ours, theirs, k)
        if off is not None:
            return off
    return None


def compare_sequences(ours: str, theirs: str) -> dict:
    """Compare a translated protein with a UniProt sequence.

    Reported statuses, from strongest to weakest agreement:
      identical                     byte-for-byte
      identical_after_met_trim      differs only by the initiator Met
      isoform_extension             one is a contiguous substring of the
                                    other: same protein, different
                                    N- or C-terminal boundary
      offset_isoform                aligned at a fixed offset with >=98%
                                    identity over the overlap
      near_identical                >=99% identity with no offset
      divergent                     genuinely different sequence

    The isoform statuses exist because UniProt canonical and GENCODE
    canonical are chosen independently and frequently disagree on the
    start codon. That is an annotation-source difference, not a
    translation error, and conflating the two hides real bugs.
    """
    base = {"len_ours": len(ours), "len_theirs": len(theirs)}
    if ours == theirs:
        return {**base, "status": "identical", "identity": 1.0,
                "n_diff": 0, "first_diff": None, "offset": 0}
    if ours.startswith("M") and ours[1:] == theirs:
        return {**base, "status": "identical_after_met_trim", "identity": 1.0,
                "n_diff": 0, "first_diff": None, "offset": 1}
    if ours and ours in theirs:
        return {**base, "status": "isoform_extension", "identity": 1.0,
                "n_diff": len(theirs) - len(ours), "first_diff": None,
                "offset": theirs.index(ours)}
    if theirs and theirs in ours:
        return {**base, "status": "isoform_extension", "identity": 1.0,
                "n_diff": len(ours) - len(theirs), "first_diff": None,
                "offset": -ours.index(theirs)}

    off = _anchor_offset(ours, theirs)
    if off is not None:
        lo = max(0, -off)
        hi = min(len(ours), len(theirs) - off)
        overlap = hi - lo
        if overlap > 0:
            diffs = [i for i in range(lo, hi) if ours[i] != theirs[i + off]]
            ident = (overlap - len(diffs)) / overlap
            if ident >= 0.98:
                return {**base, "status": "offset_isoform",
                        "identity": round(ident, 5), "n_diff": len(diffs),
                        "first_diff": (diffs[0] + 1) if diffs else None,
                        "offset": off}

    n = min(len(ours), len(theirs))
    diffs = [i for i in range(n) if ours[i] != theirs[i]]
    ident = (n - len(diffs)) / max(len(ours), len(theirs), 1)
    status = "near_identical" if ident >= 0.99 else "divergent"
    return {**base, "status": status, "identity": round(ident, 5),
            "n_diff": len(diffs) + abs(len(ours) - len(theirs)),
            "first_diff": (diffs[0] + 1) if diffs else n + 1, "offset": 0}


def best_uniprot_match(db: "UniProtDB", gene: str, seq: str):
    """The UniProt entry for `gene` that best matches `seq`.

    Picking the longest entry is wrong when a gene has several SwissProt
    entries: it can compare our protein against an unrelated paralogous
    entry and report a spurious divergence.
    """
    cands = db.for_gene(gene)
    if not cands:
        return None, None
    RANK = {"identical": 0, "identical_after_met_trim": 1,
            "isoform_extension": 2, "offset_isoform": 3,
            "near_identical": 4, "divergent": 5}
    best = best_c = None
    for e in cands:
        c = compare_sequences(seq, e.sequence)
        key = (RANK.get(c["status"], 9), -c["identity"])
        if best is None or key < best:
            best, best_c, best_e = key, c, e
    return best_e, best_c


def selenoprotein_genes(db: UniProtDB) -> set[str]:
    """Genes whose UniProt sequence contains U - these truncate at UGA
    under a naive standard-code translation and must be excluded from
    identity statistics rather than counted as failures."""
    out = set()
    for e in db.entries.values():
        if "U" in e.sequence and e.gn:
            out.update(e.gn.split())
    return out

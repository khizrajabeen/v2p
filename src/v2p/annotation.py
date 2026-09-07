"""Reference genome + GENCODE/Ensembl GTF transcript models.

Coordinate conventions used throughout this package:
  * genomic coordinates are 1-based, inclusive on both ends (GTF/VCF style)
  * transcript and CDS coordinates are 0-based offsets from the 5' end of
    the mature transcript, in transcription direction
Every conversion goes through this module so the convention lives in one
place.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

from .seqops import revcomp

_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


@dataclass
class Transcript:
    tx_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    biotype: str = ""
    tx_biotype: str = ""
    exons: list[tuple[int, int]] = field(default_factory=list)   # 1-based incl
    cds: list[tuple[int, int]] = field(default_factory=list)     # 1-based incl
    tags: set[str] = field(default_factory=set)
    cds_phase_first: int = 0   # GTF frame of the first CDS block
    sec_sites: list[int] = field(default_factory=list)  # genomic Sec codon starts

    def sort_blocks(self) -> None:
        rev = self.strand == "-"
        self.exons.sort(key=lambda x: x[0], reverse=rev)
        self.cds.sort(key=lambda x: x[0], reverse=rev)

    @property
    def is_coding(self) -> bool:
        return bool(self.cds)

    @property
    def tx_length(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)

    @property
    def start(self) -> int:
        return min(s for s, _ in self.exons)

    @property
    def end(self) -> int:
        return max(e for _, e in self.exons)

    # -- coordinate mapping ---------------------------------------------
    def genomic_to_tx(self, pos: int) -> int | None:
        """1-based genomic -> 0-based transcript offset, or None if intronic."""
        off = 0
        for s, e in self.exons:          # already in transcription order
            if s <= pos <= e:
                return off + (pos - s if self.strand == "+" else e - pos)
            off += e - s + 1
        return None

    def cds_offset_in_tx(self) -> int | None:
        """0-based transcript offset of the first translated base.

        The GTF frame column of the first CDS block is the number of
        bases to discard before the first complete codon. It is 0 for a
        normal transcript, but 1 or 2 for the ~2,000 GENCODE transcripts
        tagged `cds_start_NF`, where the 5' end of the CDS could not be
        determined. Ignoring it translates those transcripts one or two
        bases out of frame and produces a plausible-looking but entirely
        wrong protein, so the phase is applied here rather than assumed
        away.
        """
        if not self.cds:
            return None
        first = self.cds[0]
        anchor = first[0] if self.strand == "+" else first[1]
        off = self.genomic_to_tx(anchor)
        return None if off is None else off + self.cds_phase_first

    def cds_length(self) -> int:
        return sum(e - s + 1 for s, e in self.cds)

    def sec_codon_indices(self, cds_offset: int | None = None,
                          variant_tx_offset: int | None = None,
                          length_delta: int = 0) -> set[int]:
        """0-based codon indices, from the CDS start, that encode Sec.

        `variant_tx_offset` and `length_delta` shift sites that sit
        downstream of an indel. A site whose own codon is disrupted by the
        variant is dropped, since it is no longer the annotated codon.
        """
        if cds_offset is None:
            cds_offset = self.cds_offset_in_tx()
        if cds_offset is None or not self.sec_sites:
            return set()
        out: set[int] = set()
        for g in self.sec_sites:
            tx = self.genomic_to_tx(g)
            if tx is None:
                continue
            if variant_tx_offset is not None and length_delta:
                if variant_tx_offset < tx:
                    tx += length_delta
                elif tx <= variant_tx_offset < tx + 3:
                    continue          # the Sec codon itself was altered
            rel = tx - cds_offset
            if rel >= 0 and rel % 3 == 0:
                out.add(rel // 3)
        return out


class Genome:
    """Thin pyfaidx wrapper that tolerates 'chr1' vs '1' naming."""

    def __init__(self, fasta_path: str | Path):
        from pyfaidx import Fasta

        self.path = str(fasta_path)
        self.fa = Fasta(self.path, sequence_always_upper=True, as_raw=True)
        self._names = set(self.fa.keys())

    def resolve(self, chrom: str) -> str | None:
        if chrom in self._names:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in self._names:
            return alt
        for mt in ("chrM", "MT", "M", "chrMT"):
            if chrom in ("chrM", "MT", "M", "chrMT") and mt in self._names:
                return mt
        return None

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive slice on the + strand."""
        name = self.resolve(chrom)
        if name is None:
            raise KeyError(f"contig {chrom!r} not in {self.path}")
        if start > end:
            return ""
        return str(self.fa[name][start - 1:end]).upper()

    def blocks(self, chrom: str, blocks: list[tuple[int, int]], strand: str) -> str:
        """Concatenate 1-based inclusive blocks in transcription order."""
        ordered = sorted(blocks, key=lambda x: x[0])
        seq = "".join(self.fetch(chrom, s, e) for s, e in ordered)
        return revcomp(seq) if strand == "-" else seq


def _parse_attrs(field_str: str) -> dict[str, str]:
    out: dict[str, str] = {}
    tags: list[str] = []
    for k, v in _ATTR_RE.findall(field_str):
        if k == "tag":
            tags.append(v)
        else:
            out.setdefault(k, v)
    if tags:
        out["_tags"] = ",".join(tags)
    return out


class Annotation:
    """Transcript models indexed by transcript id and by gene."""

    def __init__(self) -> None:
        self.tx: dict[str, Transcript] = {}
        self.by_gene_id: dict[str, list[str]] = {}
        self.by_gene_name: dict[str, list[str]] = {}

    @classmethod
    def from_gtf(cls, gtf_path: str | Path,
                 features: tuple[str, ...] = ("exon", "CDS", "Selenocysteine"),
                 coding_only: bool = True,
                 logger=None) -> "Annotation":
        """Parse a GTF into transcript models.

        `coding_only` skips non-coding transcripts while reading. GENCODE
        v44 has ~250k transcripts but only ~90k protein-coding ones, and
        we can never emit a protein from the rest, so filtering here cuts
        both parse time and peak memory by roughly two thirds.
        """
        ann = cls()
        path = Path(gtf_path)
        opener = gzip.open if path.suffix == ".gz" else open
        n_lines = 0
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in features:
                    continue
                n_lines += 1
                a = _parse_attrs(f[8])
                tx_id = a.get("transcript_id")
                if not tx_id:
                    continue
                if coding_only:
                    tt = a.get("transcript_type", a.get("transcript_biotype", ""))
                    if tt and tt != "protein_coding":
                        continue
                t = ann.tx.get(tx_id)
                if t is None:
                    t = Transcript(
                        tx_id=tx_id,
                        gene_id=a.get("gene_id", ""),
                        gene_name=a.get("gene_name", a.get("gene_id", "")),
                        chrom=f[0],
                        strand=f[6],
                        biotype=a.get("gene_type", a.get("gene_biotype", "")),
                        tx_biotype=a.get("transcript_type",
                                         a.get("transcript_biotype", "")),
                        tags=set((a.get("_tags") or "").split(",")) - {""},
                    )
                    ann.tx[tx_id] = t
                    ann.by_gene_id.setdefault(t.gene_id, []).append(tx_id)
                    base = t.gene_id.split(".")[0]
                    if base != t.gene_id:
                        ann.by_gene_id.setdefault(base, []).append(tx_id)
                    ann.by_gene_name.setdefault(t.gene_name, []).append(tx_id)
                start, end = int(f[3]), int(f[4])
                if f[2] == "exon":
                    t.exons.append((start, end))
                elif f[2] == "Selenocysteine":
                    # GENCODE marks the exact codon that reads UGA as
                    # selenocysteine rather than as a stop.
                    t.sec_sites.append(start if t.strand == "+" else end)
                else:
                    if not t.cds:
                        t.cds_phase_first = 0 if f[7] in (".", "") else int(f[7])
                    t.cds.append((start, end))
        for t in ann.tx.values():
            t.sort_blocks()
        if logger:
            logger.info("GTF %s: %d feature lines, %d transcripts, %d genes",
                        path.name, n_lines, len(ann.tx), len(ann.by_gene_id))
        return ann

    # -- lookups ---------------------------------------------------------
    def transcripts_for_gene(self, gene: str,
                             coding_only: bool = True) -> list[Transcript]:
        ids = (self.by_gene_id.get(gene)
               or self.by_gene_id.get(gene.split(".")[0])
               or self.by_gene_name.get(gene)
               or [])
        out = [self.tx[i] for i in ids]
        if coding_only:
            out = [t for t in out if t.is_coding]
        return out

    def representative(self, gene: str,
                       coding_only: bool = True) -> Transcript | None:
        """Pick one transcript per gene, deterministically.

        Preference order: MANE_Select > Ensembl_canonical > basic tag >
        longest CDS > longest transcript > lexicographically smallest id.
        The final tie-break on the id guarantees the same choice on every
        run, which matters for reproducibility.
        """
        cands = self.transcripts_for_gene(gene, coding_only=coding_only)
        if not cands:
            return None

        def key(t: Transcript):
            return (
                0 if any(x.startswith("MANE_Select") for x in t.tags) else 1,
                0 if "Ensembl_canonical" in t.tags else 1,
                0 if "basic" in t.tags else 1,
                -t.cds_length(),
                -t.tx_length,
                t.tx_id,
            )

        return sorted(cands, key=key)[0]

    def coding_transcripts(self) -> list[Transcript]:
        return [t for t in self.tx.values() if t.is_coding]

    def ambiguous_gene_names(self) -> dict[str, list[str]]:
        """Gene names that map to more than one gene id.

        Paralogues, readthrough loci and legacy symbol collisions all
        produce these. Resolving such a name without positional context
        can pick the wrong locus, which yields either a wrong protein or,
        more often, silently no protein at all.
        """
        out: dict[str, list[str]] = {}
        for name, tids in self.by_gene_name.items():
            gids = {self.tx[t].gene_id for t in tids}
            if len(gids) > 1:
                out[name] = sorted(gids)
        return out

    def transcripts_for_gene_at(self, gene: str, chrom: str, pos: int,
                                coding_only: bool = True) -> list[Transcript]:
        """Transcripts of `gene` whose span actually contains `pos`.

        Falls back to the plain name lookup only when nothing overlaps, so
        an ambiguous symbol is resolved by locus rather than by luck.
        """
        c = chrom.lstrip("chr")
        cands = self.transcripts_for_gene(gene, coding_only=coding_only)
        hit = [t for t in cands
               if t.chrom.lstrip("chr") == c and t.start <= pos <= t.end]
        return hit or []

    def representative_at(self, gene: str, chrom: str, pos: int,
                          coding_only: bool = True) -> Transcript | None:
        """Representative transcript of `gene` at a specific locus."""
        cands = self.transcripts_for_gene_at(gene, chrom, pos,
                                             coding_only=coding_only)
        if not cands:
            return None

        def key(t: Transcript):
            return (0 if any(x.startswith("MANE_Select") for x in t.tags) else 1,
                    0 if "Ensembl_canonical" in t.tags else 1,
                    0 if "basic" in t.tags else 1,
                    -t.cds_length(), -t.tx_length, t.tx_id)

        return sorted(cands, key=key)[0]

    # -- positional lookup ----------------------------------------------
    BIN = 1 << 17          # 131 kb bins

    def build_position_index(self, logger=None) -> None:
        """Bin transcripts by genomic interval.

        A sites-only VCF carries no gene symbol, so the only way to find
        the affected transcripts is by coordinate. Binning keeps this
        O(1) per query instead of scanning ~90k transcripts 33,221 times.
        """
        self._index: dict[tuple[str, int], list[str]] = {}
        n = 0
        for t in self.tx.values():
            if not t.is_coding:
                continue
            n += 1
            key = t.chrom.lstrip("chr")
            for b in range(t.start // self.BIN, t.end // self.BIN + 1):
                self._index.setdefault((key, b), []).append(t.tx_id)
        if logger:
            logger.info("position index: %d coding transcripts in %d bins",
                        n, len(self._index))

    def transcripts_at(self, chrom: str, pos: int) -> list[Transcript]:
        """Coding transcripts whose span contains `pos`, exon-checked."""
        if not hasattr(self, "_index"):
            self.build_position_index()
        key = (chrom.lstrip("chr"), pos // self.BIN)
        out = []
        for tid in self._index.get(key, ()):
            t = self.tx[tid]
            if t.start <= pos <= t.end:
                out.append(t)
        return out

    def representatives_at(self, chrom: str, pos: int) -> list[Transcript]:
        """One transcript per gene overlapping `pos`."""
        by_gene: dict[str, list[Transcript]] = {}
        for t in self.transcripts_at(chrom, pos):
            by_gene.setdefault(t.gene_id, []).append(t)

        def key(t: Transcript):
            return (0 if any(x.startswith("MANE_Select") for x in t.tags) else 1,
                    0 if "Ensembl_canonical" in t.tags else 1,
                    0 if "basic" in t.tags else 1,
                    -t.cds_length(), -t.tx_length, t.tx_id)

        return [sorted(v, key=key)[0] for v in by_gene.values()]

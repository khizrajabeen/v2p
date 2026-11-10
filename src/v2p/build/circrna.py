"""Circular RNA: proteins translated across a back-splice junction.

A back-splice joins the 3' end of an exon to the 5' end of an exon at or
upstream of it, giving a covalently closed transcript with no free ends.
Translation can therefore run past the junction and continue around the
circle, which linear translation cannot do.

Only ORFs that cross the junction are kept. An ORF lying wholly inside
the circularised exons is also present in the linear mRNA of the same
gene, so it is not evidence of circularisation and emitting it would
claim a circRNA protein from a sequence the linear transcript already
explains. The number of crossings is recorded as `translation_cycles`:
one means the ORF reads once past the junction, more means it travelled
the circle repeatedly before meeting a stop.

The circle is laid out by repeating the sequence `max_cycles` times and
translating that, which is the same device as reading around the circle
and stopping when a stop codon appears.
"""

from __future__ import annotations

from ..annotation import Annotation, Genome, Transcript
from ..seqops import table_for_contig, translate_orf
from .smallvar import ProteinRecord

__all__ = ["CIRCRNA", "circular_sequence", "build_circrna_proteins"]

CIRCRNA = "CIRCRNA"


def circular_sequence(t: Transcript, genome: Genome, start: int,
                      end: int) -> tuple[str, list[tuple[int, int]]]:
    """Exons of `t` lying within [start, end], joined in transcription order.

    Returns (sequence, exons used). A back-splice is reported by its outer
    genomic bounds, and the circle is the exons of the host transcript
    between them, so the exon structure comes from the annotation rather
    than from the caller's coordinates.
    """
    inside = [(s, e) for s, e in t.exons if s >= start and e <= end]
    if not inside:
        return "", []
    inside.sort(key=lambda x: x[0], reverse=t.strand == "-")
    return genome.blocks(t.chrom, inside, t.strand), inside


def _orfs_crossing(seq: str, period: int, min_aa: int, require_atg: bool,
                   table: dict):
    """ORFs of a repeated circle that pass the junction at least once.

    `period` is the length of one turn, so a junction sits at every
    multiple of it. An ORF crosses when its first and last bases fall in
    different turns, and the number of turns between them is how many
    times it read past the junction.
    """
    out = []
    for frame in (0, 1, 2):
        sub = seq[frame:]
        n_full = (len(sub) // 3) * 3
        protein = translate_orf(sub[:n_full], table, to_stop=False).protein
        aa_index = 0
        for chunk in protein.split("*"):
            text, offset = chunk, aa_index
            if require_atg:
                m = text.find("M")
                text = text[m:] if m >= 0 else ""
                offset += m if m >= 0 else 0
            if text and len(text) >= min_aa:
                nt_start = frame + offset * 3
                nt_end = nt_start + len(text) * 3        # exclusive
                first_turn = nt_start // period
                last_turn = (nt_end - 1) // period
                if last_turn > first_turn:
                    out.append((text, frame, nt_start, nt_end,
                                last_turn - first_turn))
            aa_index += len(chunk) + 1                   # + the stop
    return out


def build_circrna_proteins(
    chrom: str, start: int, end: int, strand: str, gene: str,
    ann: Annotation, genome: Genome,
    min_aa: int = 30, require_atg: bool = True, max_cycles: int = 3,
    max_per_event: int = 1, source: str = "", confidence: str = "",
    logger=None,
) -> list[ProteinRecord]:
    """Proteins from one back-splice junction, longest ORF first."""
    cands = ann.transcripts_for_gene(gene, coding_only=False) if gene else []
    if not cands:
        cands = ann.transcripts_at(chrom, start) or []
    cands = [t for t in cands
             if t.chrom.lstrip("chr") == chrom.lstrip("chr")]
    if not cands:
        if logger:
            logger.warning("circRNA %s:%d-%d: no transcript for %r",
                           chrom, start, end, gene)
        return []

    out: list[ProteinRecord] = []
    seen: set[str] = set()
    for t in sorted(cands, key=lambda x: (-len(x.exons), x.tx_id)):
        circ, exons = circular_sequence(t, genome, start, end)
        if len(circ) < 3:
            continue
        table = table_for_contig(t.chrom)
        found = _orfs_crossing(circ * max_cycles, len(circ), min_aa,
                               require_atg, table)
        if not found:
            continue
        found.sort(key=lambda x: (-len(x[0]), x[1], x[2]))
        for protein, frame, nt_start, _nt_end, cycles in found[:max_per_event]:
            if protein in seen:
                continue
            seen.add(protein)
            # Where the junction falls inside this protein, 1-based: the
            # residues around it are the part no linear transcript has.
            junction_aa = max(1, (len(circ) - nt_start) // 3 + 1)
            out.append(ProteinRecord(
                seq_id=f"CIRC_{gene or t.gene_name}_{t.tx_id}_"
                       f"{chrom}_{start}_{end}",
                sequence=protein,
                variant_class=CIRCRNA,
                consequence="circrna_backsplice_orf",
                gene=gene or t.gene_name,
                transcript=t.tx_id,
                protein_change=f"p.bsj@{junction_aa}",
                variant_pos_aa=junction_aa,
                ref_protein_len=None,
                locus=f"{chrom}:{start}-{end}",
                source=source,
                confidence=str(confidence),
                novel_span=(junction_aa, len(protein)),
                notes=[f"frame{frame}", f"translation_cycles_{cycles}",
                       "crosses_backsplice_junction"]
                      + ([] if require_atg else ["no_atg_start"]),
                extra={"circ_length_nt": len(circ),
                       "n_exons": len(exons),
                       "translation_cycles": cycles,
                       "frame": frame,
                       "strand": strand or t.strand,
                       "description": f"{gene or t.gene_name} circRNA ORF "
                                      f"across the back-splice junction"},
            ))
        if out:
            break           # one host transcript is enough
    return out

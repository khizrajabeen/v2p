"""Chimeric transcript construction and translation for gene fusions.

The input fusion table gives genomic breakpoints without partner strand
or transcript context, so this module:

  1. picks a representative coding transcript for each partner gene;
  2. takes the 5' partner's mature-transcript sequence up to the
     breakpoint and the 3' partner's from the breakpoint onward - a
     breakpoint falling in an intron therefore snaps to the flanking
     exon boundary, which is the RNA-level junction that would actually
     be spliced;
  3. translates from the 5' partner's annotated start codon and reports
     whether the junction is in-frame.

Both breakpoint orientations are attempted (gene1 as 5' partner and
gene2 as 5' partner) unless the caller fixes the order, because fusion
tables do not consistently order partners 5'->3'.
"""

from __future__ import annotations

from ..annotation import Annotation, Genome, Transcript
from ..seqops import translate_orf
from .smallvar import ProteinRecord


def tx_bases_upto(t: Transcript, gpos: int) -> int:
    """Count mature-transcript bases from the 5' end through `gpos`.

    If `gpos` is intronic the count stops at the end of the last exon
    entirely 5' of it (i.e. the breakpoint snaps to the splice donor).
    """
    n = 0
    for s, e in t.exons:                     # transcription order
        if t.strand == "+":
            if e <= gpos:
                n += e - s + 1
            elif s <= gpos <= e:
                n += gpos - s + 1
                break
            else:
                break
        else:
            if s >= gpos:
                n += e - s + 1
            elif s <= gpos <= e:
                n += e - gpos + 1
                break
            else:
                break
    return n


def tx_bases_before(t: Transcript, gpos: int) -> int:
    """Mature-transcript bases strictly 5' of `gpos`."""
    n = tx_bases_upto(t, gpos)
    return n - 1 if t.genomic_to_tx(gpos) is not None else n


def _is_exonic(t: Transcript, gpos: int) -> bool:
    return t.genomic_to_tx(gpos) is not None


def build_fusion_proteins(
    gene1: str, chrom1: str, pos1: int,
    gene2: str, chrom2: str, pos2: int,
    ann: Annotation, genome: Genome,
    source: str = "", confidence: str = "",
    detection_type: str = "", validated: str = "",
    try_both_orientations: bool = True,
    min_partner_nt: int = 30,
    logger=None,
) -> list[ProteinRecord]:
    orders = [((gene1, chrom1, pos1), (gene2, chrom2, pos2))]
    if try_both_orientations:
        orders.append(((gene2, chrom2, pos2), (gene1, chrom1, pos1)))

    out: list[ProteinRecord] = []
    for (g5, c5, p5), (g3, c3, p3) in orders:
        # Resolve each partner at its own breakpoint, so an ambiguous
        # symbol cannot select a paralogue on another chromosome.
        t5 = ann.representative_at(g5, c5, p5) or ann.representative(g5)
        t3 = ann.representative_at(g3, c3, p3) or ann.representative(g3)
        if t5 is None or t3 is None:
            if logger:
                logger.warning("fusion %s--%s: no coding transcript for %s",
                               gene1, gene2, g5 if t5 is None else g3)
            continue
        if t5.chrom.lstrip("chr") != c5.lstrip("chr"):
            continue
        if t3.chrom.lstrip("chr") != c3.lstrip("chr"):
            continue

        seq5_full = genome.blocks(t5.chrom, t5.exons, t5.strand)
        seq3_full = genome.blocks(t3.chrom, t3.exons, t3.strand)

        n5 = tx_bases_upto(t5, p5)
        n3_before = tx_bases_before(t3, p3)
        head = seq5_full[:n5]
        tail = seq3_full[n3_before:]

        notes: list[str] = []
        if not _is_exonic(t5, p5):
            notes.append("bp5_intronic_snapped_to_exon_end")
        if not _is_exonic(t3, p3):
            notes.append("bp3_intronic_snapped_to_exon_start")
        if len(head) < min_partner_nt or len(tail) < min_partner_nt:
            notes.append("short_partner_segment")

        cds_off = t5.cds_offset_in_tx()
        if cds_off is None or cds_off >= len(head):
            notes.append("bp5_upstream_of_start_codon")
            continue

        chimeric = head + tail
        tr = translate_orf(chimeric[cds_off:])
        if not tr.protein:
            continue

        coding_head_nt = len(head) - cds_off
        junction_aa = coding_head_nt // 3 + 1          # 1-based residue index
        codon_intact = coding_head_nt % 3 == 0
        notes.append("junction_on_codon_boundary" if codon_intact
                     else f"junction_mid_codon(offset={coding_head_nt % 3})")

        # Does the 3' partner continue in its own native reading frame?
        # That is the definition that matters for a fusion neoantigen: the
        # downstream protein is native only if the chimeric CDS phase at the
        # junction equals the 3' partner's own CDS phase there.
        cds_off3 = t3.cds_offset_in_tx()
        native_3p: bool | None = None
        if cds_off3 is not None:
            k = n3_before - cds_off3          # 3' CDS offset of first tail base
            if k >= 0:
                native_3p = (coding_head_nt % 3) == (k % 3)
                notes.append("3p_native_frame" if native_3p
                             else "3p_frameshifted")
            else:
                notes.append("bp3_in_5UTR_of_partner")
        in_frame = bool(codon_intact and native_3p)
        if not tr.stop_found:
            notes.append("orf_runs_to_transcript_end")

        # how much of the protein is genuinely novel (3' partner derived)
        novel_start = junction_aa if in_frame else max(1, junction_aa)
        novel = (novel_start, len(tr.protein)) if len(tr.protein) >= novel_start else None

        out.append(ProteinRecord(
            seq_id=f"{g5}--{g3}_{t5.tx_id}--{t3.tx_id}_{c5}_{p5}_{c3}_{p3}",
            sequence=tr.protein,
            variant_class="FUSION",
            consequence="in_frame_fusion" if in_frame else "frameshift_fusion",
            gene=f"{g5}--{g3}",
            transcript=f"{t5.tx_id}--{t3.tx_id}",
            protein_change=f"p.junction@{junction_aa}",
            variant_pos_aa=junction_aa,
            ref_protein_len=None,
            locus=f"{c5}:{p5}::{c3}:{p3}",
            source=source,
            confidence=str(confidence),
            novel_span=novel,
            notes=notes,
            extra={"detection_type": detection_type, "validated": validated,
                   "five_prime_gene": g5, "three_prime_gene": g3,
                   "head_nt": len(head), "tail_nt": len(tail),
                   "strand5": t5.strand, "strand3": t3.strand,
                   "codon_boundary": int(codon_intact),
                   "three_prime_native_frame": (
                       "" if native_3p is None else int(native_3p))},
        ))
    return out

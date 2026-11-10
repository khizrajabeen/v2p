"""Proteoforms from alternative translation of an unchanged transcript.

Nothing here reads a variant. These forms arise from how a ribosome
reads a sequence, not from a change to it, so they are properties of the
reference annotation and are identical for every sample.

Selenocysteine termination
--------------------------
A UGA in a selenoprotein mRNA is read as selenocysteine when the SECIS
element and its machinery are engaged, and as a stop codon when they are
not. Both outcomes happen in the same cell, and the second one is the
default: incorporation is inefficient, and falls further under selenium
restriction. v2p already emits the read-through form, with `U` at the
recoded position. The terminated form is the same transcript stopping at
that codon, and is emitted here.

The truncated protein is a prefix of the read-through one, so almost all
of its tryptic peptides are shared. What is not shared is its final
peptide, which ends at the truncation instead of continuing through it.
That peptide is the evidence the terminated form was present, and it
exists in no other entry of the database.
"""

from __future__ import annotations

from ..seqops import table_for_contig, translate_orf
from .smallvar import ProteinRecord

__all__ = ["SEC_TERMINATION", "build_sec_termination_proteins"]

SEC_TERMINATION = "SEC_TERMINATION"


def build_sec_termination_proteins(transcripts, genome, min_aa: int = 8,
                                   id_prefix: str = "SECT",
                                   logger=None) -> list[ProteinRecord]:
    """One record per selenoprotein transcript, terminated at its first Sec.

    `min_aa` drops a transcript whose Sec codon sits so close to the start
    that the truncated form is too short to yield a peptide at all.
    """
    records: list[ProteinRecord] = []
    n_skipped_short = n_skipped_same = 0

    for t in sorted(transcripts, key=lambda x: x.tx_id):
        if not getattr(t, "sec_sites", None) or not t.is_coding:
            continue
        off = t.cds_offset_in_tx()
        if off is None:
            continue
        seq = genome.blocks(t.chrom, t.exons, t.strand)
        if not seq:
            continue
        table = table_for_contig(t.chrom)

        # Without the Sec annotation the UGA terminates, which is exactly
        # the form wanted here; with it, translation reads through.
        terminated = translate_orf(seq[off:], table).protein
        readthrough = translate_orf(
            seq[off:], table, sec_codons=t.sec_codon_indices(off)).protein

        if len(terminated) < min_aa:
            n_skipped_short += 1
            continue
        if terminated == readthrough:
            # No Sec codon actually recoded; nothing distinguishes them.
            n_skipped_same += 1
            continue

        gene = getattr(t, "gene_name", "") or getattr(t, "gene_id", "")
        records.append(ProteinRecord(
            seq_id=f"{id_prefix}_{gene}_{t.tx_id}",
            sequence=terminated,
            variant_class=SEC_TERMINATION,
            consequence="selenocysteine_termination",
            gene=gene,
            transcript=t.tx_id,
            protein_change=f"p.U{len(terminated) + 1}*",
            variant_pos_aa=len(terminated),
            ref_protein_len=len(readthrough),
            locus=f"{t.chrom}:{min(s for s, _ in t.exons)}"
                  f"-{max(e for _, e in t.exons)}",
            source="reference_annotation",
            novel_span=(len(terminated), len(terminated)),
            notes=["selenocysteine_termination",
                   f"readthrough_length_{len(readthrough)}"],
            extra={"n_sec_sites": len(t.sec_sites),
                   "readthrough_length": len(readthrough),
                   "description": f"{gene or t.tx_id} terminated at "
                                  f"selenocysteine codon"},
        ))

    if logger:
        logger.info("selenocysteine termination: %d record(s); skipped %d "
                    "too short and %d with no recoded codon",
                    len(records), n_skipped_short, n_skipped_same)
    return records

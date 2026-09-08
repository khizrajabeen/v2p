"""Combinatorial proteoforms: co-occurring variants applied together.

Every tool in this space, including this one until now, translates one
variant at a time. A transcript carrying three somatic SNVs yields three
proteins, each with one substitution and the other two reverted to
reference. No cell contains those proteins. The cell contains the protein
carrying all three.

That matters for a search database in a way that is easy to miss. A
tryptic peptide spanning two nearby variants exists only in the
combinatorial form: it is absent from every single-variant entry and from
the reference, so a search against a single-variant database cannot
identify it at any FDR. The closer two variants sit, the likelier they
share a peptide, and the more certain the miss.

**This is the one thing no competing tool can do**, and not because the
others are careless. Combining evidence requires ingesting it first, and
v2p is the only tool that reads DNA variants, RNA editing, fusions and
splicing into one coordinate space. A somatic SNV combined with an ADAR
edit on the same transcript is a proteoform with no purely genomic basis
*and* no purely transcriptomic basis - it needs both call sets together,
which is exactly what this pipeline already produces.

Off by default. It adds at most one entry per multi-variant transcript,
which is modest, but it makes a claim about co-occurrence that the input
does not always support - see phasing.

Phasing
-------
Co-occurrence in a call set is not phase. Two heterozygous variants may
sit on opposite alleles, in which case no single molecule carries both and
the combinatorial protein does not exist. v2p is not given phase
information, so every combined entry is a *hypothesis*, marked
`unphased` in its notes so it can be filtered. Where a caller supplies
phase, honouring it is the obvious next step; inventing it is not.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..annotation import Annotation, Genome, Transcript
from ..seqops import first_difference, revcomp, table_for_contig, translate_orf
from .smallvar import (
    ProteinRecord, _cds_offset, _is_left_anchored, _tx_interval, _tx_sequence,
)

__all__ = ["CombinableVariant", "group_by_transcript",
           "build_combinatorial_proteins", "COMBO_CLASS"]

COMBO_CLASS = "COMBO"


@dataclass(frozen=True)
class CombinableVariant:
    """One variant, reduced to what combining needs."""

    chrom: str
    pos: int
    ref: str
    alt: str
    variant_class: str = "SNV"
    source: str = ""
    confidence: str = ""

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.pos}{self.ref}>{self.alt}"


def _placement(t: Transcript, v: CombinableVariant):
    """Where `v` lands on `t`, or None.

    Returns (tx_start, tx_len, alt_allele) in transcript orientation, and
    mirrors the single-variant path including the VCF anchor rule.
    """
    e_ref, e_alt = v.ref, v.alt
    iv = _tx_interval(t, v.pos, len(v.ref))
    if iv is None and _is_left_anchored(v.ref, v.alt):
        # The anchor base is often just outside the exon while the changed
        # bases are not; same fallback as build_small_variant_proteins.
        e_ref, e_alt = v.ref[1:], v.alt[1:]
        iv = _tx_interval(t, v.pos + 1, len(e_ref)) if e_ref else None
    if iv is None:
        return None
    tx_start, tx_len = iv
    allele = e_alt if t.strand == "+" else revcomp(e_alt)
    return tx_start, tx_len, allele


def group_by_transcript(variants, ann: Annotation,
                        transcript_mode: str = "representative"):
    """transcript id -> (transcript, variants landing on it).

    Only transcripts carrying two or more variants are returned; a single
    variant is the ordinary path and is already handled.
    """
    hits: dict[str, tuple] = {}
    for v in variants:
        cands = (ann.representatives_at(v.chrom, v.pos)
                 if transcript_mode == "representative"
                 else ann.transcripts_at(v.chrom, v.pos))
        for t in cands:
            if not t.is_coding or _placement(t, v) is None:
                continue
            hits.setdefault(t.tx_id, (t, []))[1].append(v)
    return {k: val for k, val in hits.items() if len(val[1]) > 1}


def _apply_all(tx_seq: str, placed: list) -> str | None:
    """Splice every placement into the transcript.

    Applied right to left, so an earlier placement's offsets stay valid
    after a later one changes the length. Overlapping placements are
    refused rather than silently producing a mangled sequence: two
    variants claiming the same bases is a contradiction in the input, not
    a proteoform.
    """
    ordered = sorted(placed, key=lambda p: p[0], reverse=True)
    for i in range(len(ordered) - 1):
        start, _ln, _a = ordered[i]
        nxt_start, nxt_len, _b = ordered[i + 1]
        if nxt_start + nxt_len > start:
            return None
    seq = tx_seq
    for tx_start, tx_len, allele in ordered:
        seq = seq[:tx_start] + allele + seq[tx_start + tx_len:]
    return seq


def build_combinatorial_proteins(
    variants, ann: Annotation, genome: Genome,
    transcript_mode: str = "representative",
    max_variants: int = 8,
    logger=None,
) -> list[ProteinRecord]:
    """One protein per transcript carrying two or more variants.

    `max_variants` caps how many are combined on one transcript. A
    transcript with dozens of calls is far likelier to be an alignment
    artefact than a real proteoform, and translating it would put a long
    stretch of fiction into the database.
    """
    out: list[ProteinRecord] = []
    groups = group_by_transcript(variants, ann, transcript_mode)
    serial = 0

    for tx_id in sorted(groups):
        t, vs = groups[tx_id]
        if len(vs) > max_variants:
            if logger:
                logger.info("%s carries %d variants, above the cap of %d - "
                            "skipped as a likely alignment artefact",
                            tx_id, len(vs), max_variants)
            continue

        placed = [p for p in (_placement(t, v) for v in vs) if p is not None]
        if len(placed) < 2:
            continue

        tx_seq = _tx_sequence(t, genome)
        cds_off = _cds_offset(t)
        if cds_off is None:
            continue

        mut_seq = _apply_all(tx_seq, placed)
        if mut_seq is None:
            if logger:
                logger.info("%s: variants overlap, not combined", tx_id)
            continue

        # A 5'UTR indel moves the start codon, as in the single path.
        delta = sum(len(a) - ln for s, ln, a in placed if s < cds_off)
        mut_cds_off = cds_off + delta

        table = table_for_contig(t.chrom)
        ref_tr = translate_orf(tx_seq[cds_off:], table,
                               sec_codons=t.sec_codon_indices(cds_off))
        alt_tr = translate_orf(mut_seq[mut_cds_off:], table)
        if not alt_tr.protein or alt_tr.protein == ref_tr.protein:
            continue

        idx = first_difference(ref_tr.protein, alt_tr.protein)
        vs_sorted = sorted(vs, key=lambda v: (v.chrom, v.pos))
        classes = sorted({v.variant_class for v in vs_sorted})
        serial += 1

        notes = ["combinatorial", "unphased"]
        if len(classes) > 1:
            # The case no other tool can reach: evidence types combined.
            notes.append("cross_evidence_" + "+".join(classes))

        out.append(ProteinRecord(
            seq_id=f"COMBO_{t.tx_id}_{serial:06d}",
            sequence=alt_tr.protein,
            variant_class=COMBO_CLASS,
            consequence="combinatorial",
            gene=t.gene_name or t.gene_id,
            transcript=t.tx_id,
            protein_change=f"p.[{len(vs_sorted)}_changes]",
            variant_pos_aa=idx + 1,
            ref_protein_len=len(ref_tr.protein),
            locus=";".join(v.locus for v in vs_sorted),
            source=";".join(sorted({v.source for v in vs_sorted if v.source})),
            novel_span=(idx + 1, len(alt_tr.protein)),
            notes=notes,
            extra={"n_variants": len(vs_sorted),
                   "loci": [v.locus for v in vs_sorted],
                   "classes": classes,
                   "cross_evidence": len(classes) > 1,
                   "description": f"{t.gene_name or t.tx_id} combinatorial "
                                  f"proteoform ({len(vs_sorted)} co-occurring "
                                  f"variants)"},
        ))
    return out

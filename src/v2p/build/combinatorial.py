"""Combinatorial proteoforms: co-occurring variants applied together.

This pipeline, until now, translated one variant at a time. A transcript
carrying three somatic SNVs yielded three proteins, each with one
substitution and the other two reverted to reference. No cell contains
those proteins. The cell contains the protein carrying all three.

Why it matters, and how much
----------------------------
A tryptic peptide spanning two nearby variants exists only in the
combinatorial form: it is absent from every single-variant entry and from
the reference, so a search against a single-variant database cannot
identify it at any FDR. This is not a hypothetical corner. ProHap
(Nature Methods 2024, doi:10.1038/s41592-024-02506-0) measured it for
common germline haplotypes and found that **12.42% of discoverable amino
acid substitutions can share a tryptic peptide with another
substitution**. Those peptides are unsearchable in a single-variant
database.

What is and is not new here
---------------------------
Combining is not new, and this module did not invent it: ProHap's whole
design is combination, from *phased* genotypes, and its companion ProVar
takes sample-level VCFs but "considers each allele independently". This
module works in ProVar's scope - one sample, calls that are often
unphased - and combines there.

The part that is genuinely unavailable elsewhere is combination *across
evidence types*. A somatic SNV co-occurring with an A-to-I edit on one
transcript cannot be represented by any haplotype panel: RNA editing is
not in the genome and never appears in a VCF of genotypes. That is a
limit of the input a genotype-based tool is given, not an oversight in
its design.

Every entry must earn its place
-------------------------------
Database size is not free. Variant-peptide FDR is already the weak point
of proteogenomics, and every added sequence costs discriminating power,
so an entry contributing no new peptide is worse than useless. A
combinatorial protein is therefore emitted only when at least one of its
tryptic peptides is absent from the reference protein *and* from every
single-variant protein of the same transcript. Two variants a thousand
residues apart produce a protein whose every peptide is already in the
database; that entry is dropped, and the number dropped is logged.

Phasing
-------
Co-occurrence in a call set is not phase. Two heterozygous variants may
sit on opposite alleles, in which case no molecule carries both and the
combinatorial protein does not exist.

Where the caller reports phase (`GT` with `|` and a `PS` phase set) it is
honoured: variants combine within a haplotype, entries are marked
`phased`, and two variants on opposite haplotypes are never combined -
that combination is not uncertain, it is false. A homozygous-alt variant
sits on both haplotypes and joins both groups.

Where phase is absent or partial the combination is a *hypothesis*,
marked `unphased` so it can be filtered downstream. Protein changes use
HGVS allele notation, which distinguishes the two cases: `p.[A12V;G45S]`
when phased, `p.[A12V(;)G45S]` when not.

Off by default.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..annotation import Annotation, Genome, Transcript
from ..peptides import digest
from ..seqops import first_difference, revcomp, table_for_contig, translate_orf
from .smallvar import (
    ProteinRecord, _cds_offset, _is_left_anchored, _tx_interval, _tx_sequence,
)

__all__ = ["CombinableVariant", "group_by_transcript", "phase_groups",
           "cooccurring_peptides", "build_combinatorial_proteins",
           "COMBO_CLASS"]

COMBO_CLASS = "COMBO"


@dataclass(frozen=True)
class CombinableVariant:
    """One variant, reduced to what combining needs.

    `phase_set` is the caller's PS tag and `haplotypes` the haplotype
    indices carrying the ALT allele - {0}, {1}, or {0, 1} when
    homozygous. Both empty means the caller reported no phase.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    variant_class: str = "SNV"
    source: str = ""
    confidence: str = ""
    phase_set: str = ""
    haplotypes: frozenset = field(default_factory=frozenset)

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.pos}{self.ref}>{self.alt}"

    @property
    def is_phased(self) -> bool:
        return bool(self.phase_set) and bool(self.haplotypes)


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


def phase_groups(variants) -> list[tuple[str, bool, list]]:
    """Split one transcript's variants into sets that can co-exist.

    Returns (label, phased, members) tuples, keeping only those with two
    or more members. A phased group is a fact about one haplotype, and
    each is emitted separately rather than merged, because two haplotypes
    are two molecules.

    An unphased group is a hypothesis, and is formed only when the
    caller has not already settled the question:

    - some variant has no GT, or no PS -> nothing is known about how it
      sits relative to the others
    - the variants span more than one PS block -> each block is phased
      within itself, but their relative phase is unknown

    It is deliberately *not* formed when every variant sits in one PS
    block, because there the caller has said which haplotype each is on.
    Two variants it placed on opposite haplotypes are not uncertain:
    no molecule carries both, and combining them would invent one.
    """
    groups: list[tuple[str, bool, list]] = []
    by_hap: dict[tuple[str, int], list] = {}
    blocks: set[str] = set()
    any_unphased = False
    for v in variants:
        if v.is_phased:
            blocks.add(v.phase_set)
            for h in sorted(v.haplotypes):
                by_hap.setdefault((v.phase_set, h), []).append(v)
        else:
            any_unphased = True
    for (ps, h), members in sorted(by_hap.items()):
        groups.append((f"{ps}|{h}", True, members))
    if any_unphased or len(blocks) > 1 or not by_hap:
        groups.append(("unphased", False, list(variants)))
    return [g for g in groups if len(g[2]) > 1]


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


def _translate(tx_seq: str, placed: list, cds_off: int, t: Transcript):
    """Apply `placed` to the transcript and translate. None if refused."""
    mut = _apply_all(tx_seq, placed)
    if mut is None:
        return None
    # A 5'UTR indel moves the start codon, as in the single-variant path.
    delta = sum(len(a) - ln for s, ln, a in placed if s < cds_off)
    tr = translate_orf(mut[cds_off + delta:], table_for_contig(t.chrom))
    return tr.protein or None


def _aa_change(ref_p: str, alt_p: str) -> str:
    """The change between two proteins, in short HGVS-like form."""
    if not alt_p or ref_p == alt_p:
        return "="
    i = first_difference(ref_p, alt_p)
    if i is None or i >= len(ref_p):
        return "ext*?"
    was = ref_p[i]
    if len(ref_p) == len(alt_p) and i < len(alt_p):
        return f"{was}{i + 1}{alt_p[i]}"
    if i >= len(alt_p):
        return f"{was}{i + 1}*"
    return f"{was}{i + 1}fs"


def _allele_changes(ref_p: str, alt_p: str, limit: int = 8) -> list[str]:
    """Residue differences between the reference and the combined form.

    Read off the combined protein rather than from each variant in turn,
    because two variants can share a codon. Describing those separately
    would claim two different residues at one position - `D225N` and
    `D225E` - when the protein carries a third residue that is neither.
    Variants that are synonymous, alone or together, contribute no
    residue and so do not appear.
    """
    if len(ref_p) != len(alt_p):
        return [_aa_change(ref_p, alt_p)]
    diffs = [f"{ref_p[i]}{i + 1}{alt_p[i]}"
             for i in range(len(ref_p)) if ref_p[i] != alt_p[i]]
    if len(diffs) > limit:
        return diffs[:limit] + [f"+{len(diffs) - limit}_more"]
    return diffs or ["="]


def cooccurring_peptides(combo: str, singles, ref: str,
                         missed_cleavages: int = 2, min_len: int = 6,
                         max_len: int = 50,
                         equate_il: bool = True) -> set[str]:
    """Peptides of `combo` that no single-variant form and no reference has.

    These exist only because the variants co-occur. If this set is empty
    the combinatorial entry contributes nothing a search could find, and
    adding it inflates the database - which costs FDR power - for no
    gain.

    Isoleucine and leucine are folded together because they are isobaric:
    a peptide differing only by I/L is not detectably new.
    """
    def _key(p: str) -> str:
        return p.replace("I", "L") if equate_il else p

    covered = set()
    for s in [ref, *singles]:
        if s:
            covered |= {_key(p) for p in
                        digest(s, missed_cleavages, min_len, max_len)}
    return {p for p in digest(combo, missed_cleavages, min_len, max_len)
            if _key(p) not in covered}


def build_combinatorial_proteins(
    variants, ann: Annotation, genome: Genome,
    transcript_mode: str = "representative",
    max_variants: int = 8,
    missed_cleavages: int = 2,
    require_new_peptide: bool = True,
    allow_unphased: bool = True,
    logger=None,
) -> list[ProteinRecord]:
    """One protein per haplotype carrying two or more variants.

    `max_variants` caps how many are combined on one transcript. A
    transcript with dozens of calls is far likelier to be an alignment
    artefact than a real proteoform, and translating it would put a long
    stretch of fiction into the database.

    `require_new_peptide` drops any combination whose every tryptic
    peptide already exists in a single-variant entry or the reference.
    Turn it off only to inspect what is being dropped.

    `allow_unphased` admits combinations the caller did not phase. It is
    on by default because sites-only and unphased VCFs are the common
    case - HCC1395's truth set has no sample column at all - but every
    entry records which it is, so the hypotheses can be filtered out
    later without rebuilding the database.
    """
    out: list[ProteinRecord] = []
    groups = group_by_transcript(variants, ann, transcript_mode)
    serial = 0
    n_no_peptide = n_overlap = n_capped = n_unphased_skipped = 0

    for tx_id in sorted(groups):
        t, tx_vars = groups[tx_id]
        cds_off = _cds_offset(t)
        if cds_off is None:
            continue
        tx_seq = _tx_sequence(t, genome)
        ref_protein = translate_orf(
            tx_seq[cds_off:], table_for_contig(t.chrom),
            sec_codons=t.sec_codon_indices(cds_off)).protein
        seen_seq: set[str] = set()

        for label, phased, vs in phase_groups(tx_vars):
            if not phased and not allow_unphased:
                n_unphased_skipped += 1
                continue
            if len(vs) > max_variants:
                n_capped += 1
                if logger:
                    logger.info("%s %s carries %d variants, above the cap of "
                                "%d - skipped as a likely alignment artefact",
                                tx_id, label, len(vs), max_variants)
                continue

            vs_sorted = sorted(vs, key=lambda v: (v.chrom, v.pos))
            pairs = [(v, _placement(t, v)) for v in vs_sorted]
            pairs = [(v, p) for v, p in pairs if p is not None]
            if len(pairs) < 2:
                continue

            combo = _translate(tx_seq, [p for _v, p in pairs], cds_off, t)
            if combo is None:
                n_overlap += 1
                if logger:
                    logger.info("%s %s: variants overlap, not combined",
                                tx_id, label)
                continue
            if combo == ref_protein or combo in seen_seq:
                continue

            singles = [_translate(tx_seq, [p], cds_off, t) or ""
                       for _v, p in pairs]
            changes = _allele_changes(ref_protein, combo)

            new_peps = cooccurring_peptides(combo, singles, ref_protein,
                                            missed_cleavages)
            if require_new_peptide and not new_peps:
                n_no_peptide += 1
                continue

            seen_seq.add(combo)
            serial += 1
            classes = sorted({v.variant_class for v, _p in pairs})
            sep = ";" if phased else "(;)"
            idx = first_difference(ref_protein, combo) or 0

            notes = ["combinatorial", "phased" if phased else "unphased"]
            if len(classes) > 1:
                # The case no other tool can reach: evidence types combined.
                notes.append("cross_evidence_" + "+".join(classes))

            out.append(ProteinRecord(
                seq_id=f"COMBO_{t.tx_id}_{serial:06d}",
                sequence=combo,
                variant_class=COMBO_CLASS,
                consequence="combinatorial",
                gene=t.gene_name or t.gene_id,
                transcript=t.tx_id,
                protein_change=f"p.[{sep.join(changes)}]",
                variant_pos_aa=idx + 1,
                ref_protein_len=len(ref_protein),
                locus=";".join(v.locus for v, _p in pairs),
                source=";".join(sorted({v.source for v, _p in pairs
                                        if v.source})),
                novel_span=(idx + 1, len(combo)),
                notes=notes,
                extra={"n_variants": len(pairs),
                       "loci": [v.locus for v, _p in pairs],
                       "classes": classes,
                       "cross_evidence": len(classes) > 1,
                       "phased": phased,
                       "phase": "phased" if phased else "unphased",
                       "phase_set": label if phased else "",
                       "n_cooccurring_peptides": len(new_peps),
                       "cooccurring_peptides": sorted(new_peps)[:20],
                       "description": (
                           f"{t.gene_name or t.tx_id} combinatorial "
                           f"proteoform ({len(pairs)} co-occurring "
                           f"variants, "
                           f"{'phased' if phased else 'unphased'})")},
            ))

    if logger and (n_no_peptide or n_overlap or n_capped
                   or n_unphased_skipped):
        logger.info("combinatorial: dropped %d contributing no new peptide, "
                    "%d overlapping, %d over the variant cap, %d unphased "
                    "(--allow-unphased is off)",
                    n_no_peptide, n_overlap, n_capped, n_unphased_skipped)
    return out

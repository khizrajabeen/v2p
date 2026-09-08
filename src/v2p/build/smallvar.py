"""Apply a small variant (SNV, MNV, InDel, RNA-editing site) to a
transcript and re-translate.

Design note
-----------
The variant is applied to the *mature transcript* sequence rather than to
genomic sequence. That is what makes one code path serve both DNA
variants and RNA-editing sites: an A-to-I edit is, at the transcript
level, indistinguishable from a genomic A>G substitution, and a
minus-strand gene's genomic T>C becomes the same transcript-level A>G
automatically once the exon blocks are reverse-complemented.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..annotation import Annotation, Genome, Transcript
from ..nmd import predict_nmd
from ..seqops import (
    Translation, first_difference, revcomp, table_for_contig, translate_orf,
)


@dataclass
class ProteinRecord:
    """One variant protein sequence plus everything needed for its header."""

    seq_id: str
    sequence: str
    variant_class: str
    consequence: str
    gene: str
    transcript: str
    protein_change: str = ""
    variant_pos_aa: int | None = None       # 1-based residue of first change
    ref_protein_len: int | None = None
    locus: str = ""
    source: str = ""
    confidence: str = ""
    novel_span: tuple[int, int] | None = None   # 1-based inclusive, novel AA
    notes: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)


# Many variants land on the same transcript - 33k variants across ~20k
# genes means heavy repetition - and rebuilding a mature transcript costs
# one pyfaidx fetch per exon. Cache by transcript id; the cache is keyed
# on the Genome object too, so switching reference cannot serve stale
# sequence.
_TX_SEQ_CACHE: dict[tuple[int, str], str] = {}
_TX_OFF_CACHE: dict[tuple[int, str], int | None] = {}
_CACHE_MAX = 60000


def clear_sequence_cache() -> None:
    _TX_SEQ_CACHE.clear()
    _TX_OFF_CACHE.clear()


def _tx_sequence(t: Transcript, genome: Genome) -> str:
    key = (id(genome), t.tx_id)
    seq = _TX_SEQ_CACHE.get(key)
    if seq is None:
        seq = genome.blocks(t.chrom, t.exons, t.strand)
        if len(_TX_SEQ_CACHE) >= _CACHE_MAX:
            _TX_SEQ_CACHE.clear()
        _TX_SEQ_CACHE[key] = seq
    return seq


def _cds_offset(t: Transcript) -> int | None:
    key = (0, t.tx_id)
    if key not in _TX_OFF_CACHE:
        if len(_TX_OFF_CACHE) >= _CACHE_MAX:
            _TX_OFF_CACHE.clear()
        _TX_OFF_CACHE[key] = t.cds_offset_in_tx()
    return _TX_OFF_CACHE[key]


def _tx_interval(t: Transcript, gpos: int, glen: int) -> tuple[int, int] | None:
    """Transcript offsets covered by genomic [gpos, gpos+glen-1].

    Returns (tx_start, tx_len) or None when the interval is not fully
    exonic / not contiguous on the transcript (i.e. it crosses an intron).
    """
    a = t.genomic_to_tx(gpos)
    b = t.genomic_to_tx(gpos + glen - 1)
    if a is None or b is None:
        return None
    lo, hi = (a, b) if a <= b else (b, a)
    if hi - lo + 1 != glen:
        return None
    return lo, glen


def classify_consequence(ref_p: str, alt_p: str, ref_nt_len: int,
                         alt_nt_len: int, alt_stop_found: bool) -> str:
    d = alt_nt_len - ref_nt_len
    if d % 3 != 0:
        return "frameshift"
    if ref_p == alt_p:
        return "synonymous"
    if alt_p and ref_p and alt_p[0] != ref_p[0] and ref_p[0] == "M" and alt_p[0] != "M":
        return "start_lost"
    if d > 0:
        return "inframe_insertion"
    if d < 0:
        return "inframe_deletion"
    if len(alt_p) < len(ref_p) and alt_stop_found:
        return "stop_gained"
    if len(alt_p) > len(ref_p):
        return "stop_lost"
    return "missense"


def build_small_variant_proteins(
    chrom: str, pos: int, ref: str, alt: str,
    ann: Annotation, genome: Genome,
    variant_class: str = "SNV",
    genes: str = "",
    source: str = "",
    confidence: str = "",
    transcript_mode: str = "representative",
    logger=None,
) -> list[ProteinRecord]:
    """Return one ProteinRecord per affected coding transcript.

    `transcript_mode`:
      'representative' - one transcript per gene (MANE/canonical preferred)
      'all'            - every coding transcript overlapping the variant
    """
    out: list[ProteinRecord] = []

    # Candidate transcripts. Gene symbols come free with an ANNOVAR table
    # but a sites-only VCF has none, so fall back to a coordinate lookup.
    # Positional lookup is also the safer answer where a symbol is present
    # but stale, so it is used whenever the symbol resolves to nothing.
    cands: list[Transcript] = []
    gene_list = [x.strip() for x in genes.replace(",", ";").split(";")
                 if x.strip() and x.strip() not in (".", "NONE", "NA")]
    for g in gene_list:
        # Resolve the symbol at this locus. 59 UniProt symbols and a
        # comparable number of GENCODE gene names are ambiguous; taking
        # the first match can land on a paralogue elsewhere, which then
        # produces no protein at all rather than an obvious error.
        if transcript_mode == "representative":
            t = ann.representative_at(g, chrom, pos) or ann.representative(g)
            if t:
                cands.append(t)
        else:
            at = ann.transcripts_for_gene_at(g, chrom, pos)
            cands.extend(at or ann.transcripts_for_gene(g))

    if not cands:
        span_end = pos + max(len(ref), 1) - 1
        if transcript_mode == "representative":
            cands = ann.representatives_at(chrom, pos)
            if not cands and span_end != pos:
                cands = ann.representatives_at(chrom, span_end)
        else:
            cands = ann.transcripts_at(chrom, pos)
            if not cands and span_end != pos:
                cands = ann.transcripts_at(chrom, span_end)

    # de-duplicate, keep deterministic order
    seen: set[str] = set()
    cands = [t for t in cands if not (t.tx_id in seen or seen.add(t.tx_id))]
    cands.sort(key=lambda t: t.tx_id)

    for t in cands:
        if genome.resolve(t.chrom) is None or genome.resolve(chrom) is None:
            continue
        if t.chrom.lstrip("chr") != chrom.lstrip("chr"):
            continue
        if not t.is_coding:
            continue

        iv = _tx_interval(t, pos, len(ref))
        if iv is None:
            continue                      # intronic, or spans a splice site
        tx_start, tx_len = iv

        tx_seq = _tx_sequence(t, genome)
        cds_off = _cds_offset(t)
        if cds_off is None:
            continue

        obs_ref = tx_seq[tx_start:tx_start + tx_len]
        exp_ref = ref if t.strand == "+" else revcomp(ref)
        note: list[str] = []
        if obs_ref != exp_ref:
            note.append(f"REF_MISMATCH(expected={exp_ref},found={obs_ref})")
            if logger:
                logger.warning(
                    "REF mismatch %s:%d %s>%s on %s: transcript has %s",
                    chrom, pos, ref, alt, t.tx_id, obs_ref)

        alt_tx_allele = alt if t.strand == "+" else revcomp(alt)
        mut_seq = tx_seq[:tx_start] + alt_tx_allele + tx_seq[tx_start + tx_len:]

        # a 5'UTR indel shifts where the CDS begins
        delta = len(alt_tx_allele) - tx_len
        mut_cds_off = cds_off + (delta if tx_start < cds_off else 0)
        if tx_start < cds_off and delta:
            note.append("5UTR_indel_shifts_CDS_start")

        # chrM uses the vertebrate mitochondrial code; GENCODE annotates
        # which UGA codons are selenocysteine rather than stop.
        codon_table = table_for_contig(t.chrom)
        ref_sec = t.sec_codon_indices(cds_off)
        alt_sec = t.sec_codon_indices(cds_off, variant_tx_offset=tx_start,
                                      length_delta=delta)
        ref_tr: Translation = translate_orf(tx_seq[cds_off:], codon_table,
                                            sec_codons=ref_sec)
        alt_tr: Translation = translate_orf(mut_seq[mut_cds_off:], codon_table,
                                            sec_codons=alt_sec)
        if ref_sec:
            note.append(f"selenoprotein_{len(ref_sec)}_Sec_codons")
        if t.chrom in ("chrM", "chrMT", "MT", "M"):
            note.append("mitochondrial_genetic_code")

        if not alt_tr.protein:
            note.append("empty_alt_protein")

        # A UTR variant is exonic, so it reaches this point, but it lies
        # outside the translated region and therefore leaves the protein
        # unchanged. Without this check it is indistinguishable from a
        # synonymous coding change, which inverts the nonsynonymous-to-
        # synonymous ratio and makes the recovery accounting wrong.
        orf_nt = 3 * len(ref_tr.protein) + (3 if ref_tr.stop_found else 0)
        cds_end_tx = cds_off + orf_nt
        var_last = tx_start + tx_len - 1
        if var_last < cds_off:
            cons = "5_prime_UTR"
        elif tx_start >= cds_end_tx:
            cons = "3_prime_UTR"
        else:
            cons = classify_consequence(
                ref_tr.protein, alt_tr.protein,
                ref_nt_len=tx_len, alt_nt_len=len(alt_tx_allele),
                alt_stop_found=alt_tr.stop_found,
            )
        if not ref_tr.stop_found:
            note.append("ref_orf_no_stop")
        if not alt_tr.stop_found:
            note.append("alt_orf_runs_to_transcript_end")
        if ref_tr.ambiguous_codons or alt_tr.ambiguous_codons:
            note.append("ambiguous_bases_translated_as_X")

        idx = first_difference(ref_tr.protein, alt_tr.protein)
        pos_aa = idx + 1 if cons != "synonymous" else None
        pchange = ""
        if cons == "synonymous":
            pchange = "p.(=)"
        elif cons == "frameshift" and idx < len(ref_tr.protein):
            new_aa = alt_tr.protein[idx] if idx < len(alt_tr.protein) else "*"
            ext = len(alt_tr.protein) - idx
            pchange = f"p.{ref_tr.protein[idx]}{idx+1}{new_aa}fs*{ext}"
        elif cons == "stop_gained" and idx < len(ref_tr.protein):
            # A premature stop truncates the protein, so the alt sequence
            # simply ends here and `idx` lands past its end. Without this
            # branch it fell through to the `del` case below and produced
            # `p.Q70del` - a single-residue deletion - for a nonsense
            # variant, contradicting its own CSQ=stop_gained and
            # misreporting the consequence to anyone reading the header.
            # HGVS spells this `p.Gln70Ter`, one-letter `p.Q70*`.
            pchange = f"p.{ref_tr.protein[idx]}{idx+1}*"
        elif idx < len(ref_tr.protein) and idx < len(alt_tr.protein):
            pchange = f"p.{ref_tr.protein[idx]}{idx+1}{alt_tr.protein[idx]}"
        elif idx < len(ref_tr.protein):
            pchange = f"p.{ref_tr.protein[idx]}{idx+1}del"

        novel = None
        if cons == "frameshift":
            novel = (idx + 1, len(alt_tr.protein))
        elif pos_aa:
            novel = (pos_aa, pos_aa)

        # Nonsense-mediated decay: a truncating variant whose stop sits
        # well upstream of the last junction is probably degraded, so the
        # protein is unlikely to exist. Flag rather than drop.
        nmd = {"nmd": "not_applicable"}
        if cons in ("frameshift", "stop_gained"):
            ex_lengths = [e - s_ + 1 for s_, e in t.exons]
            stop_off = mut_cds_off + 3 * len(alt_tr.protein)
            nmd = predict_nmd(
                stop_codon_tx_offset=stop_off,
                exon_lengths=ex_lengths,
                variant_tx_offset=tx_start,
                length_delta=delta,
                stop_found=alt_tr.stop_found,
            )
            if nmd["nmd"] == "likely":
                note.append("NMD_likely")

        out.append(ProteinRecord(
            seq_id=f"{t.gene_name}_{t.tx_id}_{chrom}_{pos}_{ref}_{alt}",
            sequence=alt_tr.protein,
            variant_class=variant_class,
            consequence=cons,
            gene=t.gene_name,
            transcript=t.tx_id,
            protein_change=pchange,
            variant_pos_aa=pos_aa,
            ref_protein_len=len(ref_tr.protein),
            locus=f"{chrom}:{pos}{ref}>{alt}",
            source=source,
            confidence=str(confidence),
            novel_span=novel,
            notes=note,
            extra={"strand": t.strand, "tx_biotype": t.tx_biotype,
                   "ref_protein": ref_tr.protein,
                   "nmd": nmd.get("nmd"), "nmd_reason": nmd.get("reason", ""),
                   "nmd_distance": nmd.get("distance_to_last_junction")},
        ))
    return out

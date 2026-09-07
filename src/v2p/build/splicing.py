"""Alternative-splicing events -> isoform protein sequences.

Two strategies, applied in order:

  A. Junction matching (preferred). Each SUPPA2 event defines the intron
     junctions that distinguish the two alternative forms. Annotated
     coding transcripts of the gene that contain the required junctions
     are selected and translated as-is. Nothing is invented: the
     sequences are real GENCODE isoforms, which keeps false-positive
     peptides out of the search space.

  B. De novo construction (fallback, SE / RI / A3 / A5 only). If no
     annotated transcript carries the junction - which is exactly the
     interesting case for a tumour-specific event - the exon chain of the
     representative transcript is edited to realise the event and then
     translated. These records are flagged `constructed=1` so they can be
     filtered separately downstream.

SUPPA2 event id coordinate grammar implemented here:
  SE  <e1>-<s2> : <e2>-<s3>
  A5  <e1>-<s3> : <e2>-<s3>
  A3  <e1>-<s2> : <e1>-<s3>
  RI  <s1> : <e1>-<s2> : <e2>
  MX  <e1>-<s2> : <e2>-<s4> : <e1>-<s3> : <e3>-<s4>
  AF  <s1> : <e1>-<s3> : <s2> : <e2>-<s3>
  AL  <e1>-<s2> : <e2> : <e1>-<s3> : <e3>
Coordinates are genomic and 1-based regardless of strand.
"""

from __future__ import annotations

from typing import Any

from ..annotation import Annotation, Genome, Transcript
from ..seqops import first_atg_offset, table_for_contig, translate_orf
from .smallvar import ProteinRecord

Junction = tuple[int, int]   # (donor_end, acceptor_start), genomic, 1-based


def transcript_junctions(t: Transcript) -> set[Junction]:
    """Intron junctions as (last base of exon, first base of next exon)."""
    blocks = sorted(t.exons, key=lambda x: x[0])
    return {(blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1)}


def event_junctions(event_type: str,
                    groups: list[list[int]]) -> dict[str, list[Junction]] | None:
    """Required junctions for each alternative form of an event.

    Keys are 'form1' (SUPPA's first/inclusion form) and 'form2'.
    Returns None if the group arity does not match the event grammar,
    which guards against silently mis-reading a malformed event id.
    """
    g = groups
    try:
        if event_type == "SE" and len(g) == 2:
            e1, s2 = g[0]
            e2, s3 = g[1]
            return {"form1": [(e1, s2), (e2, s3)], "form2": [(e1, s3)]}
        if event_type == "A5" and len(g) == 2:
            e1, s3 = g[0]
            e2, _ = g[1]
            return {"form1": [(e1, s3)], "form2": [(e2, s3)]}
        if event_type == "A3" and len(g) == 2:
            e1, s2 = g[0]
            _, s3 = g[1]
            return {"form1": [(e1, s2)], "form2": [(e1, s3)]}
        if event_type == "RI" and len(g) == 3:
            e1, s2 = g[1]
            return {"form1": [], "form2": [(e1, s2)]}
        if event_type == "MX" and len(g) == 4:
            e1, s2 = g[0]
            e2, s4 = g[1]
            _, s3 = g[2]
            e3, _ = g[3]
            return {"form1": [(e1, s2), (e2, s4)],
                    "form2": [(e1, s3), (e3, s4)]}
        if event_type in ("AF", "AL") and len(g) == 4:
            # SUPPA2 emits AF/AL coordinate groups in opposite order on the
            # two strands: [single, pair, single, pair] one way and
            # [pair, single, pair, single] the other. Positional unpacking
            # therefore drops every event on one strand - 234 of 1027 in
            # the HCC1395 set. Select the two *pair* groups instead, which
            # are the distinguishing junctions regardless of orientation.
            pairs = [x for x in g if len(x) == 2]
            if len(pairs) != 2:
                return None
            (a1, a2), (b1, b2) = pairs
            return {"form1": [(min(a1, a2), max(a1, a2))],
                    "form2": [(min(b1, b2), max(b1, b2))]}
    except (ValueError, IndexError):
        return None
    return None


def _translate_transcript(t: Transcript, genome: Genome) -> tuple[str, list[str]]:
    seq = genome.blocks(t.chrom, t.exons, t.strand)
    off = t.cds_offset_in_tx()
    notes: list[str] = []
    if off is None:
        off = first_atg_offset(seq)
        if off is None:
            return "", ["no_orf_found"]
        notes.append("start_codon_inferred_first_ATG")
    tr = translate_orf(seq[off:], table_for_contig(t.chrom),
                       sec_codons=t.sec_codon_indices(off))
    if not tr.stop_found:
        notes.append("orf_runs_to_transcript_end")
    if tr.ambiguous_codons:
        notes.append("ambiguous_bases_translated_as_X")
    if tr.selenocysteines:
        notes.append(f"selenoprotein_{tr.selenocysteines}_Sec")
    return tr.protein, notes


# --------------------------------------------------------------------------
# Strategy B: de novo isoform construction
# --------------------------------------------------------------------------

def _construct_exons(event_type: str, groups: list[list[int]],
                     template: Transcript) -> list[tuple[int, int]] | None:
    """Edit the template exon chain to realise the *included* form."""
    ex = sorted(template.exons, key=lambda x: x[0])
    g = groups
    if event_type == "SE" and len(g) == 2:
        e1, s2 = g[0]
        e2, s3 = g[1]
        # insert cassette exon (s2, e2) between the flanking exons
        new = [b for b in ex if b[1] <= e1 or b[0] >= s3]
        new.append((s2, e2))
        return sorted(set(new), key=lambda x: x[0])
    if event_type == "RI" and len(g) == 3:
        s1 = g[0][0]
        e1, s2 = g[1]
        e2 = g[2][0]
        new = [b for b in ex if b[1] < s1 or b[0] > e2]
        new.append((s1, e2))                      # intron retained
        return sorted(set(new), key=lambda x: x[0])
    if event_type == "A3" and len(g) == 2:
        e1, s2 = g[0]
        _, s3 = g[1]
        new = []
        for s, e in ex:
            if s in (s2, s3):
                new.append((min(s2, s3), e))      # use the distal acceptor
            else:
                new.append((s, e))
        return sorted(set(new), key=lambda x: x[0])
    if event_type == "A5" and len(g) == 2:
        e1, _ = g[0]
        e2, _ = g[1]
        new = []
        for s, e in ex:
            if e in (e1, e2):
                new.append((s, max(e1, e2)))      # use the distal donor
            else:
                new.append((s, e))
        return sorted(set(new), key=lambda x: x[0])
    return None


def build_splicing_proteins(
    event_id: str, event_type: str, gene_id: str, chrom: str, strand: str,
    coord_groups: list[list[int]],
    ann: Annotation, genome: Genome,
    source: str = "", confidence: str = "",
    dpsi: float | None = None, p_val: float | None = None,
    allow_construction: bool = True,
    max_tx_per_form: int = 2,
    logger=None,
) -> list[ProteinRecord]:
    out: list[ProteinRecord] = []
    reqs = event_junctions(event_type, coord_groups)
    if reqs is None:
        if logger:
            logger.warning("event %s: %s grammar has %d coordinate groups, "
                           "cannot interpret", event_id, event_type,
                           len(coord_groups))
        return out

    cands = ann.transcripts_for_gene(gene_id, coding_only=True)
    if not cands:
        if logger:
            logger.warning("event %s: no coding transcript for gene %s",
                           event_id, gene_id)
        return out

    jx_by_tx = {t.tx_id: transcript_junctions(t) for t in cands}
    matched_any = False

    for form, required in reqs.items():
        hits = [t for t in cands if set(required) <= jx_by_tx[t.tx_id]]
        # RI form1 (intron retained) has no junction to require - detect it
        # positively instead, as an exon spanning the whole intron
        if event_type == "RI" and form == "form1":
            e1, s2 = coord_groups[1]
            hits = [t for t in cands
                    if any(s <= e1 and e >= s2 for s, e in t.exons)]
        if not hits:
            continue
        matched_any = True

        def rank(t: Transcript):
            return (0 if any(x.startswith("MANE_Select") for x in t.tags) else 1,
                    0 if "Ensembl_canonical" in t.tags else 1,
                    -t.cds_length(), t.tx_id)

        for t in sorted(hits, key=rank)[:max_tx_per_form]:
            prot, notes = _translate_transcript(t, genome)
            if not prot:
                continue
            out.append(ProteinRecord(
                seq_id=f"{gene_id}_{event_type}_{form}_{t.tx_id}",
                sequence=prot,
                variant_class=f"AS_{event_type}",
                consequence=f"{event_type}_{form}_annotated_isoform",
                gene=t.gene_name or gene_id,
                transcript=t.tx_id,
                protein_change="",
                variant_pos_aa=None,
                ref_protein_len=len(prot),
                locus=f"{chrom}:{coord_groups[0][0]}",
                source=source,
                confidence=str(confidence),
                novel_span=None,
                notes=notes + [f"form={form}", "matched_annotated_transcript"],
                extra={"event_id": event_id, "event_type": event_type,
                       "strand": strand, "dPSI": dpsi, "p_val": p_val,
                       "constructed": 0, "gtf_chrom": t.chrom,
                       "gtf_strand": t.strand, "gtf_exons": list(t.exons),
                       "gtf_cds": list(t.cds)},
            ))

    if matched_any or not allow_construction:
        return out

    # ---- fallback: build the isoform ourselves --------------------------
    template = ann.representative(gene_id)
    if template is None:
        return out
    exons = _construct_exons(event_type, coord_groups, template)
    if not exons:
        if logger:
            logger.info("event %s (%s): no annotated match and no "
                        "construction rule", event_id, event_type)
        return out

    synth = Transcript(
        tx_id=f"{template.tx_id}_{event_type}synth", gene_id=template.gene_id,
        gene_name=template.gene_name, chrom=template.chrom,
        strand=template.strand, biotype=template.biotype,
        tx_biotype="constructed_isoform", exons=list(exons), cds=[],
    )
    synth.sort_blocks()
    seq = genome.blocks(synth.chrom, synth.exons, synth.strand)
    off = first_atg_offset(seq)
    if off is None:
        return out
    tr = translate_orf(seq[off:])
    if not tr.protein:
        return out
    out.append(ProteinRecord(
        seq_id=f"{gene_id}_{event_type}_constructed_{template.tx_id}",
        sequence=tr.protein,
        variant_class=f"AS_{event_type}",
        consequence=f"{event_type}_constructed_isoform",
        gene=template.gene_name or gene_id,
        transcript=synth.tx_id,
        protein_change="",
        variant_pos_aa=None,
        ref_protein_len=None,
        locus=f"{chrom}:{coord_groups[0][0]}",
        source=source,
        confidence=str(confidence),
        novel_span=None,
        notes=["constructed_isoform", "start_codon_inferred_first_ATG"]
              + ([] if tr.stop_found else ["orf_runs_to_transcript_end"]),
        extra={"event_id": event_id, "event_type": event_type,
               "strand": strand, "dPSI": dpsi, "p_val": p_val,
               "constructed": 1, "template_tx": template.tx_id,
               "gtf_chrom": synth.chrom, "gtf_strand": synth.strand,
               "gtf_exons": list(synth.exons),
               "gtf_cds": [], "cds_tx_offset": off},
    ))
    return out

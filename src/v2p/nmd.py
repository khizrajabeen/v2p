"""Nonsense-mediated decay prediction.

A premature termination codon more than ~50 nt upstream of the final
exon-exon junction is recognised by the EJC-dependent surveillance
machinery and the transcript is degraded, so the truncated protein is
usually never made in appreciable amounts. Emitting those sequences
without a flag inflates a search database with proteoforms that do not
exist; removing them outright would discard the NMD-escaping minority
that genuinely do. So they are emitted and labelled.

The 50-nt boundary is the standard heuristic, not a law. Single-exon
transcripts have no junction and always escape. Stops in the last exon
escape. Long 3' UTRs and a handful of characterised transcripts escape
for other reasons this does not model, which is why the output says
`likely` rather than `yes`.
"""

from __future__ import annotations

NMD_BOUNDARY_NT = 50


def junction_offsets(exon_lengths: list[int]) -> list[int]:
    """Transcript offsets of each exon-exon junction, 0-based.

    A junction offset is the transcript coordinate of the first base of
    the following exon. A single-exon transcript yields an empty list.
    """
    out: list[int] = []
    total = 0
    for ln in exon_lengths[:-1]:
        total += ln
        out.append(total)
    return out


def predict_nmd(stop_codon_tx_offset: int | None,
                exon_lengths: list[int],
                variant_tx_offset: int | None = None,
                length_delta: int = 0,
                stop_found: bool = True,
                boundary: int = NMD_BOUNDARY_NT) -> dict:
    """Classify a termination codon as NMD-triggering or NMD-escaping.

    `stop_codon_tx_offset` is the 0-based transcript offset of the first
    base of the stop codon **in the variant transcript**. An indel
    upstream of a junction shifts every downstream junction by
    `length_delta`, so junctions are corrected before comparison.
    """
    if stop_codon_tx_offset is None or not stop_found:
        return {"nmd": "no_stop", "distance_to_last_junction": None,
                "n_junctions": max(0, len(exon_lengths) - 1)}

    junctions = junction_offsets(exon_lengths)
    if length_delta and variant_tx_offset is not None:
        junctions = [j + length_delta if j > variant_tx_offset else j
                     for j in junctions]

    if not junctions:
        return {"nmd": "escapes", "reason": "single_exon_transcript",
                "distance_to_last_junction": None, "n_junctions": 0}

    last = junctions[-1]
    distance = last - stop_codon_tx_offset
    if distance > boundary:
        return {"nmd": "likely", "reason": f"stop_{distance}nt_upstream_of_"
                f"last_junction", "distance_to_last_junction": distance,
                "n_junctions": len(junctions)}
    if stop_codon_tx_offset >= last:
        reason = "stop_in_last_exon"
    else:
        reason = f"stop_within_{boundary}nt_of_last_junction"
    return {"nmd": "escapes", "reason": reason,
            "distance_to_last_junction": distance,
            "n_junctions": len(junctions)}

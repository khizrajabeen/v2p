"""Non-canonical ORFs: three-frame translation of non-coding transcripts.

lncRNAs, pseudogenes and untranslated regions are not annotated as coding,
but ribosome profiling and proteomics keep finding peptides from them. This
is the one capability where pypgatk is genuinely ahead of v2p, so it exists
for parity.

**These entries are excluded by default and must stay that way.** Three
frames across every non-coding transcript multiplies database size several
fold, and an inflated search space costs sensitivity: more candidate
peptides means a higher score threshold at the same FDR, which loses real
identifications elsewhere. `--include-noncanonical` is an informed choice,
not a default.

Nothing here reads a variant. These ORFs are properties of the reference
annotation, so they are identical for every sample; they are emitted
alongside variant records rather than derived from them.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..seqops import CODON_TABLE
from .smallvar import ProteinRecord

__all__ = [
    "NC_CLASSES", "NoncanonicalORF", "classify_biotype",
    "find_orfs", "build_noncanonical_proteins", "build_utr_orfs",
    "drop_known_proteins",
    "DEFAULT_MIN_AA",
]

# The default floor. 30 residues is roughly the shortest stretch that can
# yield an identifying tryptic peptide at all; below it an entry adds
# search space without adding evidence.
DEFAULT_MIN_AA = 30

NC_LNCRNA = "NC_LNCRNA"
NC_PSEUDOGENE = "NC_PSEUDOGENE"
NC_UTR = "NC_UTR"
NC_ALTORF = "NC_ALTORF"

NC_CLASSES = (NC_LNCRNA, NC_PSEUDOGENE, NC_UTR, NC_ALTORF)

# GENCODE biotype -> our class. Matched on substring because the vocabulary
# is long and keeps growing (`lncRNA`, `processed_pseudogene`,
# `transcribed_unitary_pseudogene`, ...); an unrecognised non-coding
# biotype falls back to NC_ALTORF rather than being dropped silently.
_BIOTYPE_HINTS = (
    ("pseudogene", NC_PSEUDOGENE),
    ("lncrna", NC_LNCRNA),
    ("lincrna", NC_LNCRNA),
    ("antisense", NC_LNCRNA),
    ("processed_transcript", NC_LNCRNA),
    ("retained_intron", NC_LNCRNA),
)

# Non-coding but not protein candidates: structural or regulatory RNAs,
# far too short and far too numerous to be worth three-frame translating.
_SKIP_BIOTYPES = ("mirna", "snrna", "snorna", "rrna", "trna", "misc_rna",
                  "scarna", "srna", "vault_rna", "ribozyme")


def classify_biotype(biotype: str) -> str | None:
    """Map a transcript biotype to an NC class, or None to skip it."""
    b = (biotype or "").strip().lower()
    if not b:
        return None
    if any(s in b for s in _SKIP_BIOTYPES):
        return None
    if b == "protein_coding":
        return None
    for needle, cls in _BIOTYPE_HINTS:
        if needle in b:
            return cls
    return NC_ALTORF


@dataclass(frozen=True)
class NoncanonicalORF:
    """One open reading frame found in a transcript."""

    protein: str
    frame: int             # 0, 1 or 2, on the transcript's own strand
    nt_start: int          # 0-based offset into the mature transcript
    nt_end: int            # exclusive
    started_at_atg: bool
    stop_found: bool

    @property
    def length(self) -> int:
        return len(self.protein)


def _translate_frame(nt: str, frame: int):
    """Yield (protein, aa_offset_in_frame, stop_found) for each ORF."""
    seq = nt[frame:]
    n = (len(seq) // 3) * 3
    aas = [CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, n, 3)]
    run: list[str] = []
    start = 0
    for i, aa in enumerate(aas):
        if aa == "*":
            if run:
                yield "".join(run), start, True
            run = []
            start = i + 1
            continue
        if not run:
            start = i
        run.append(aa)
    if run:
        # No terminating stop: the ORF runs off the end of the transcript.
        yield "".join(run), start, False


def find_orfs(nt: str, min_aa: int = DEFAULT_MIN_AA,
              require_atg: bool = True) -> list[NoncanonicalORF]:
    """Three-frame ORFs on one strand, longest first.

    `require_atg` trims each ORF to its first methionine, which is the
    conservative reading: an ORF reported from an arbitrary start codon is
    a claim about translation initiation that the sequence alone cannot
    support. With it off, the full inter-stop stretch is kept and the frame
    is recorded so the caller can judge.
    """
    out: list[NoncanonicalORF] = []
    for frame in (0, 1, 2):
        for prot, aa_off, stop_found in _translate_frame(nt, frame):
            started_at_atg = False
            if require_atg:
                m = prot.find("M")
                if m < 0:
                    continue
                aa_off += m
                prot = prot[m:]
                started_at_atg = True
            if len(prot) < min_aa:
                continue
            nt_start = frame + aa_off * 3
            out.append(NoncanonicalORF(
                protein=prot, frame=frame, nt_start=nt_start,
                nt_end=nt_start + len(prot) * 3,
                started_at_atg=started_at_atg, stop_found=stop_found))
    # Longest first, then frame and position, so the order is stable across
    # runs - the reproducibility test depends on it.
    out.sort(key=lambda o: (-o.length, o.frame, o.nt_start))
    return out


def drop_known_proteins(records, reference_sequences, logger=None):
    """Remove ORFs whose protein is already in the reference proteome.

    A three-frame translation of a pseudogene routinely reproduces the
    protein of the parent gene, and an ORF that duplicates a sequence the
    database already contains is pure search-space inflation: it cannot
    yield a peptide the reference does not already explain, but it does
    raise the score threshold at a fixed FDR for everything else.

    Matching is on the exact sequence. A one-residue difference is a real
    proteoform and is kept.
    """
    known = {s for s in reference_sequences if s}
    kept, dropped = [], 0
    for r in records:
        if r.sequence in known:
            dropped += 1
            continue
        kept.append(r)
    if logger and dropped:
        logger.info("dropped %d non-canonical ORF(s) identical to a "
                    "reference protein", dropped)
    return kept, dropped


def build_utr_orfs(
    transcripts, genome, min_aa: int = DEFAULT_MIN_AA,
    require_atg: bool = True, max_per_transcript: int = 1,
    id_prefix: str = "NC",
) -> list[ProteinRecord]:
    """ORFs in the UTRs of *coding* transcripts - uORFs and dORFs.

    Distinct from `build_noncanonical_proteins`, which handles whole
    non-coding transcripts. A 5-prime uORF is the better studied case:
    they regulate translation of the main ORF and their peptides do turn
    up in MS data, but they are invisible to any pipeline that only
    translates the annotated CDS.

    The main CDS itself is excluded - it is already emitted by the normal
    path, and re-translating it here would duplicate every reference
    protein under a second class.
    """
    records: list[ProteinRecord] = []
    serial = 0
    for tx in sorted(transcripts, key=lambda t: t.tx_id):
        if not getattr(tx, "cds", None) or not tx.exons:
            continue                     # non-coding has no UTR to speak of
        start = tx.cds_offset_in_tx()
        if start is None:
            continue
        cds_len = sum(e - s + 1 for s, e in tx.cds)
        nt = genome.blocks(tx.chrom, tx.exons, tx.strand)
        if not nt:
            continue

        regions = (("5_prime", nt[:start]), ("3_prime", nt[start + cds_len:]))
        for which, sub in regions:
            if len(sub) < min_aa * 3:
                continue
            for orf in find_orfs(sub, min_aa=min_aa,
                                 require_atg=require_atg)[:max_per_transcript]:
                serial += 1
                gene = (getattr(tx, "gene_name", "")
                        or getattr(tx, "gene_id", ""))
                notes = [f"frame{orf.frame}", f"{which}_utr_orf",
                         "three_frame_translation"]
                if not orf.stop_found:
                    notes.append("orf_runs_to_region_end")
                records.append(ProteinRecord(
                    seq_id=f"{id_prefix}_{NC_UTR}_{serial:06d}",
                    sequence=orf.protein,
                    variant_class=NC_UTR,
                    consequence="noncanonical_orf",
                    gene=gene,
                    transcript=tx.tx_id,
                    locus=f"{tx.chrom}:{min(s for s, _ in tx.exons)}"
                          f"-{max(e for _, e in tx.exons)}",
                    source="reference_annotation",
                    novel_span=(1, orf.length),
                    notes=notes,
                    extra={"frame": orf.frame, "utr": which,
                           "biotype": getattr(tx, "tx_biotype", ""),
                           "description": f"{gene or tx.tx_id} {which} UTR "
                                          f"ORF (frame {orf.frame})"},
                ))
    return records


def build_noncanonical_proteins(
    transcripts, genome, min_aa: int = DEFAULT_MIN_AA,
    require_atg: bool = True, max_per_transcript: int = 1,
    id_prefix: str = "NC",
) -> list[ProteinRecord]:
    """Three-frame translate non-coding transcripts into ProteinRecords.

    `max_per_transcript` caps how many ORFs each transcript contributes.
    The default of 1 keeps the longest, because emitting every ORF in every
    frame is what makes this feature dangerous to database size.
    """
    records: list[ProteinRecord] = []
    serial = 0
    # Sorted so output order does not depend on dict iteration order.
    for tx in sorted(transcripts, key=lambda t: t.tx_id):
        if getattr(tx, "cds", None):
            continue                      # it is coding; not our business
        cls = classify_biotype(getattr(tx, "tx_biotype", "")
                               or getattr(tx, "biotype", ""))
        if cls is None or not tx.exons:
            continue

        nt = genome.blocks(tx.chrom, tx.exons, tx.strand)
        if not nt:
            continue

        for orf in find_orfs(nt, min_aa=min_aa,
                             require_atg=require_atg)[:max_per_transcript]:
            serial += 1
            gene = getattr(tx, "gene_name", "") or getattr(tx, "gene_id", "")
            notes = [f"frame{orf.frame}", "three_frame_translation"]
            if not orf.stop_found:
                notes.append("orf_runs_to_transcript_end")
            if not orf.started_at_atg:
                notes.append("no_atg_start")
            records.append(ProteinRecord(
                seq_id=f"{id_prefix}_{cls}_{serial:06d}",
                sequence=orf.protein,
                variant_class=cls,
                consequence="noncanonical_orf",
                gene=gene,
                transcript=tx.tx_id,
                protein_change="",
                locus=f"{tx.chrom}:{min(s for s, _ in tx.exons)}"
                      f"-{max(e for _, e in tx.exons)}",
                source="reference_annotation",
                confidence="",
                novel_span=(1, orf.length),
                notes=notes,
                extra={"frame": orf.frame,
                       "biotype": getattr(tx, "tx_biotype", ""),
                       "orf_nt_start": orf.nt_start,
                       "orf_nt_end": orf.nt_end,
                       "description": f"{gene or tx.tx_id} non-canonical ORF "
                                      f"(frame {orf.frame})"},
            ))
    return records

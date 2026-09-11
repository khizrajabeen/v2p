#!/usr/bin/env python3
"""How often do two substitutions land in one tryptic peptide?

ProHap (Nature Methods 2024, doi:10.1038/s41592-024-02506-0) reports
that for **12.42%** of discoverable amino acid substitutions encoded by
common germline haplotypes, two or more co-occur in the same peptide
after tryptic digestion. Those peptides exist in no single-variant
database, which is the entire argument for combining.

That figure is for common germline haplotypes, where one transcript
routinely carries several common variants. Somatic calls are far
sparser, so the same measurement on a tumour sample should come out
lower. This script measures it, whatever it turns out to be.

Definitions, kept to the ProHap wording:

- a substitution is **discoverable** if it falls inside at least one
  tryptic peptide of detectable length
- it **co-occurs** if at least one such peptide also contains another
  substitution

The digest runs on the protein carrying *all* of the transcript's
substitutions, not on the reference, because that is the molecule whose
peptides a search engine would have to match - and because a
substitution can create or destroy a cleavage site itself (HCC1395's
FKTN p.[D225K] introduces a lysine, and with it a new cleavage site).

Usage:
  python3 benchmarks/cooccurrence_stat.py \
      --manifest results/tables/unified_variant_manifest.tsv \
      --genome ref/GRCh38.primary_assembly.genome.fa \
      --gtf ref/gencode.v44.annotation.gtf.gz \
      --out benchmarks/cooccurrence_hcc1395.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.annotation import Annotation, Genome            # noqa: E402
from v2p.build.combinatorial import (                    # noqa: E402
    CombinableVariant, _apply_all, _placement,
)
from v2p.build.smallvar import (                         # noqa: E402
    _cds_offset, _tx_sequence,
)
from v2p.seqops import table_for_contig, translate_orf   # noqa: E402

SMALL_CLASSES = {"SNV", "MNV", "INDEL", "RNA_EDITING"}
_TRYPSIN = re.compile(r"(?<=[KR])(?!P)")


def peptide_spans(seq: str, missed_cleavages: int = 2,
                  min_len: int = 6, max_len: int = 50):
    """Tryptic peptides as (start, end) residue offsets, end exclusive.

    The same rule as v2p.peptides.digest - cleave after K or R, not
    before P - but keeping positions, which is what deciding "are these
    two substitutions in one peptide" needs.
    """
    seq = seq.replace("*", "")
    if not seq:
        return []
    cuts = [0] + [m.start() for m in _TRYPSIN.finditer(seq)] + [len(seq)]
    spans = []
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + missed_cleavages + 2, len(cuts))):
            ln = cuts[j] - cuts[i]
            if min_len <= ln <= max_len:
                spans.append((cuts[i], cuts[j]))
            elif ln > max_len:
                break
    return spans


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--transcript-mode", default="representative",
                    choices=["all", "representative"])
    ap.add_argument("--missed-cleavages", type=int, default=2)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    ann = Annotation.from_gtf(Path(args.gtf))
    genome = Genome(Path(args.genome))

    variants: list[CombinableVariant] = []
    with open(args.manifest, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["variant_class"] not in SMALL_CLASSES:
                continue
            p = json.loads(row["payload_json"])
            # Substitutions only. An indel shifts every downstream
            # residue, so "which peptide is it in" is a different
            # question, and the ProHap figure is about substitutions.
            if len(p["ref"]) != 1 or len(p["alt"]) != 1:
                continue
            variants.append(CombinableVariant(
                p["chrom"], p["pos"], p["ref"], p["alt"],
                row["variant_class"], row["source"]))

    by_tx: dict[str, tuple] = {}
    for v in variants:
        cands = (ann.representatives_at(v.chrom, v.pos)
                 if args.transcript_mode == "representative"
                 else ann.transcripts_at(v.chrom, v.pos))
        for t in cands:
            if not t.is_coding or _placement(t, v) is None:
                continue
            by_tx.setdefault(t.tx_id, (t, []))[1].append(v)

    rows = []
    n_sub = n_disc = n_cooc = 0
    for tx_id in sorted(by_tx):
        t, vs = by_tx[tx_id]
        cds_off = _cds_offset(t)
        if cds_off is None:
            continue
        tx_seq = _tx_sequence(t, genome)
        placed = [(v, _placement(t, v)) for v in vs]
        placed = [(v, p) for v, p in placed if p is not None]
        mut = _apply_all(tx_seq, [p for _v, p in placed])
        if mut is None:                       # overlapping calls
            continue
        prot = translate_orf(mut[cds_off:],
                             table_for_contig(t.chrom)).protein
        if not prot:
            continue

        # Residue index of each substitution. Substitutions do not
        # change length, so transcript offsets map straight to codons.
        res = set()
        for _v, (tx_start, _ln, _a) in placed:
            if tx_start < cds_off:
                continue                      # 5'UTR, not translated
            idx = (tx_start - cds_off) // 3
            if idx < len(prot):
                res.add(idx)
        if not res:
            continue

        spans = peptide_spans(prot, args.missed_cleavages)
        disc, cooc = set(), set()
        for start, end in spans:
            inside = {i for i in res if start <= i < end}
            disc |= inside
            if len(inside) > 1:
                cooc |= inside

        n_sub += len(res)
        n_disc += len(disc)
        n_cooc += len(cooc)
        rows.append({"transcript": tx_id, "gene": t.gene_name or t.gene_id,
                     "n_subs": len(res), "n_discoverable": len(disc),
                     "n_cooccurring": len(cooc)})

    pct = (100.0 * n_cooc / n_disc) if n_disc else 0.0
    print(f"transcripts with >=1 coding substitution : {len(rows)}")
    print(f"substitutions on those transcripts       : {n_sub}")
    print(f"discoverable (inside a tryptic peptide)  : {n_disc}")
    print(f"of those, sharing a peptide with another : {n_cooc}")
    print(f"co-occurrence rate                       : {pct:.2f}%")
    print("ProHap, common germline haplotypes       : 12.42%")

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t",
                               fieldnames=["transcript", "gene", "n_subs",
                                           "n_discoverable", "n_cooccurring"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"per-transcript detail -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

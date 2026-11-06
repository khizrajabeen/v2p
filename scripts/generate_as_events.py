#!/usr/bin/env python3
"""Generate SUPPA2-style alternative splicing event ids from an annotation.

    python scripts/generate_as_events.py --gtf ref/gencode.v44.annotation.gtf.gz \
        --out events/gencode_v44_events.tsv

Writes one row per event: the id, its type, gene, strand, and the junctions
each of its two forms requires. The ids follow SUPPA2's grammar, which is
what v2p's splicing path reads.

Coordinates come from v2p's own Annotation, so the junctions written here
and the junctions v2p derives from the ids are produced by one convention
rather than two that have to agree.

Event grammar, coordinates genomic and 1-based regardless of strand:
    SE  <e1>-<s2> : <e2>-<s3>
    A5  <e1>-<s3> : <e2>-<s3>      pairs sharing the acceptor
    A3  <e1>-<s2> : <e1>-<s3>      pairs sharing the donor
    RI  <s1> : <e1>-<s2> : <e2>
    MX  <e1>-<s2> : <e2>-<s4> : <e1>-<s3> : <e3>-<s4>
    AF  <s1> : <e1>-<s3> : <s2> : <e2>-<s3>
    AL  <e1>-<s2> : <e2> : <e1>-<s3> : <e3>

A5 and A3 are the same shape and are told apart by which side of the
junction the two forms disagree on, which depends on the strand: the donor
is the lower coordinate on '+' and the higher one on '-'.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.annotation import Annotation, Transcript        # noqa: E402

Junction = tuple[int, int]      # (lower, upper), genomic, 1-based


def junctions(t: Transcript) -> list[Junction]:
    """Introns of one transcript as (last exon base, next exon base)."""
    ex = sorted(t.exons, key=lambda x: x[0])
    return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]


def _fmt(groups: list[str], etype: str, chrom: str, strand: str,
         gene: str) -> str:
    return f"{gene};{etype}:{chrom}:" + ":".join(groups) + f":{strand}"


def _se(by_tx, all_jx):
    """Cassette exon: both flanking junctions in one form, the span in another."""
    out = []
    for tx_id, jx in by_tx.items():
        for i in range(len(jx) - 1):
            e1, s2 = jx[i]
            e2, s3 = jx[i + 1]
            if (e1, s3) in all_jx:              # something skips the exon
                out.append((f"{e1}-{s2}", f"{e2}-{s3}"))
    return out


def _alt_ss(all_jx, strand):
    """Junction pairs sharing one end: alternative donor (A5) or acceptor (A3)."""
    by_lo, by_hi = defaultdict(set), defaultdict(set)
    for lo, hi in all_jx:
        by_lo[lo].add(hi)
        by_hi[hi].add(lo)
    a5, a3 = [], []
    # The two forms are written shortest intron first, on both strands and
    # for both types. Measured against 31,912 SUPPA ids: sorting the
    # variable coordinate numerically happens to give that order when the
    # shared coordinate is the lower one, and the reverse of it when the
    # shared coordinate is the upper one.
    for lo, his in by_lo.items():                  # shared lower coordinate
        for a, b in _pairs(sorted(his)):           # ascending = shorter first
            groups = (f"{lo}-{a}", f"{lo}-{b}")
            # On '+' the upper coordinate is the acceptor -> A3.
            (a3 if strand == "+" else a5).append(groups)
    for hi, los in by_hi.items():                  # shared upper coordinate
        for a, b in _pairs(sorted(los, reverse=True)):   # descending
            groups = (f"{a}-{hi}", f"{b}-{hi}")
            # On '+' the lower coordinate is the donor -> A5.
            (a5 if strand == "+" else a3).append(groups)
    return a5, a3


def _pairs(xs):
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            yield xs[i], xs[j]


def _ri(tx_exons, all_jx):
    """Retained intron: an exon spanning a junction another transcript splices."""
    out = []
    for exons in tx_exons.values():
        for s, e in exons:
            for j1, j2 in all_jx:
                if s < j1 and e > j2:
                    out.append((f"{s}", f"{j1}-{j2}", f"{e}"))
    return out


def _mx(by_tx, all_jx):
    """Mutually exclusive exons: two cassettes that never co-occur."""
    out = []
    # Collect (upstream donor, exon, downstream acceptor) triples.
    triples = defaultdict(set)
    for jx in by_tx.values():
        for i in range(len(jx) - 1):
            e1, s2 = jx[i]
            e2, s3 = jx[i + 1]
            triples[(e1, s3)].add((s2, e2))
    for (e1, s4), exons in triples.items():
        if len(exons) < 2:
            continue
        for (s2, e2), (s3, e3) in _pairs(sorted(exons)):
            if e2 >= s3:                       # overlapping: not exclusive
                continue
            out.append((f"{e1}-{s2}", f"{e2}-{s4}", f"{e1}-{s3}", f"{e3}-{s4}"))
    return out


def _first_last(transcripts, strand):
    """Alternative first and last exons.

    Each form is the terminal exon's outer coordinate plus the junction
    tying it to the exon both forms share. SUPPA writes the single and
    the paired group in opposite order on the two strands, and orders the
    two forms by the variable junction coordinate, ascending.
    """
    af, al = [], []
    firsts, lasts = set(), set()
    for t in transcripts:
        ex = sorted(t.exons, key=lambda x: x[0])
        if len(ex) < 2:
            continue
        if strand == "+":
            # first exon is leftmost: donor moves, acceptor is shared
            firsts.add((ex[0][1], ex[0][0], ex[1][0]))     # var, outer, shared
            lasts.add((ex[-1][0], ex[-1][1], ex[-2][1]))   # var, outer, shared
        else:
            # transcription runs right to left
            firsts.add((ex[-1][0], ex[-1][1], ex[-2][1]))
            lasts.add((ex[0][1], ex[0][0], ex[1][0]))

    for (v1, o1, k1), (v2, o2, k2) in _pairs(sorted(firsts)):
        if k1 != k2 or v1 == v2:
            continue
        if strand == "+":
            af.append((f"{o1}", f"{v1}-{k1}", f"{o2}", f"{v2}-{k1}"))
        else:
            af.append((f"{k1}-{v1}", f"{o1}", f"{k1}-{v2}", f"{o2}"))
    for (v1, o1, k1), (v2, o2, k2) in _pairs(sorted(lasts)):
        if k1 != k2 or v1 == v2:
            continue
        if strand == "+":
            al.append((f"{k1}-{v1}", f"{o1}", f"{k1}-{v2}", f"{o2}"))
        else:
            al.append((f"{o1}", f"{v1}-{k1}", f"{o2}", f"{v2}-{k1}"))
    return af, al


def generate(ann: Annotation, coding_only: bool = False, logger=None):
    """Yield (event_id, etype, gene_id, chrom, strand, n_transcripts)."""
    by_gene: dict[str, list[Transcript]] = defaultdict(list)
    for t in ann.tx.values():
        if coding_only and not t.is_coding:
            continue
        by_gene[t.gene_id].append(t)

    seen: set[str] = set()
    for gene in sorted(by_gene):
        txs = by_gene[gene]
        if len(txs) < 2:
            continue
        chrom, strand = txs[0].chrom, txs[0].strand
        by_tx = {t.tx_id: junctions(t) for t in txs}
        tx_exons = {t.tx_id: sorted(t.exons, key=lambda x: x[0]) for t in txs}
        all_jx = {j for jx in by_tx.values() for j in jx}
        if not all_jx:
            continue

        emit = []
        emit += [("SE", g) for g in _se(by_tx, all_jx)]
        a5, a3 = _alt_ss(all_jx, strand)
        emit += [("A5", g) for g in a5]
        emit += [("A3", g) for g in a3]
        emit += [("RI", g) for g in _ri(tx_exons, all_jx)]
        emit += [("MX", g) for g in _mx(by_tx, all_jx)]
        af, al = _first_last(txs, strand)
        emit += [("AF", g) for g in af]
        emit += [("AL", g) for g in al]

        for etype, groups in emit:
            eid = _fmt(list(groups), etype, chrom, strand, gene)
            if eid in seen:
                continue
            seen.add(eid)
            yield eid, etype, gene, chrom, strand, len(txs)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--coding-only", action="store_true",
                    help="only transcripts with a CDS; v2p translates no "
                         "others, so events confined to non-coding "
                         "transcripts can never yield a protein")
    ap.add_argument("--types", default="SE,A3,A5,RI,MX,AF,AL")
    args = ap.parse_args(argv)

    keep = {t.strip() for t in args.types.split(",") if t.strip()}
    print("loading annotation ...", flush=True)
    ann = Annotation.from_gtf(args.gtf, coding_only=args.coding_only)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = defaultdict(int)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["event_id", "event_type", "gene_id", "chrom", "strand",
                    "n_transcripts"])
        for row in generate(ann, coding_only=args.coding_only):
            if row[1] not in keep:
                continue
            counts[row[1]] += 1
            w.writerow(row)

    total = sum(counts.values())
    print(f"{total:,} events -> {args.out}")
    for t in sorted(counts):
        print(f"  {t}  {counts[t]:>9,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

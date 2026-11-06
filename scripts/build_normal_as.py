#!/usr/bin/env python3
"""Select the AS events of a catalogue that are observed in normal tissue.

    python scripts/build_normal_as.py \
        --events events/gencode_v44_events.tsv \
        --gtex   incoming/normal_junctions/gtex/GTEx_..._junctions.gct.gz \
        --recount3 incoming/normal_junctions/recount3 \
        --outdir data/normal_as

Writes one CSV per source, in the shape v2p's splicing path reads.

An event is kept when every junction it requires is observed. The
junctions come from v2p's own `event_junctions`, so the test applied here
and the matching v2p does later are the same function, not two readings
of the same grammar.

Only junctions some event needs are tracked, which keeps one pass over
each source and a few million keys in memory instead of tens of millions.

All three sources state introns: start is the first intron base, end the
last. v2p states the flanking exon bases, so a junction converts as
(start - 1, end + 1). Verified against GENCODE v44 for Snaptron, the
GTEx .gct and the recount3 catalogues alike.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.build.splicing import event_junctions              # noqa: E402
from v2p.parse.inputs import parse_suppa_event_id           # noqa: E402

Junction = tuple[int, int]


def required_junctions(event_id: str) -> tuple[str, set[tuple[str, int, int]]]:
    """(event_type, junctions the event needs) keyed by chromosome."""
    parsed = parse_suppa_event_id(event_id)
    if parsed is None:
        return "", set()
    reqs = event_junctions(parsed["event_type"], parsed["coord_groups"])
    if reqs is None:
        return parsed["event_type"], set()
    chrom = parsed["chrom"]
    out = {(chrom, a, b) for js in reqs.values() for a, b in js}
    return parsed["event_type"], out


def load_events(path: Path):
    """Yield (row, event_type, required junctions)."""
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            etype, jx = required_junctions(row["event_id"])
            yield row, etype, jx


def scan_gtex(gct: Path, wanted: set) -> dict:
    """Junctions of the GTEx .gct, counted once each.

    Only the row name is parsed. The matrix is 357,746 junctions wide by
    17,382 samples; summing every cell would mean six billion values for
    a question the row's existence already answers, since GTEx publishes
    an analysis-grade junction set rather than every aligner artefact.
    """
    hits: dict = {}
    n = 0
    with gzip.open(gct, "rb") as fh:
        fh.readline(), fh.readline(), fh.readline()      # gct header
        for raw in fh:
            cut = raw.find(b"\t")
            if cut < 0:
                continue
            n += 1
            name = raw[:cut].decode()
            chrom, _, rest = name.partition("_")
            start, _, end = rest.partition("_")
            try:
                key = (chrom, int(start) - 1, int(end) + 1)
            except ValueError:
                continue
            if key in wanted:
                hits[key] = hits.get(key, 0) + 1
    print(f"  GTEx: {n:,} junctions scanned, {len(hits):,} of the "
          f"{len(wanted):,} wanted were found", flush=True)
    return hits


def scan_recount3(folder: Path, wanted: set, pattern: str) -> dict:
    """Count, per junction, how many tissue catalogues list it."""
    hits: dict = defaultdict(int)
    files = sorted(folder.glob(pattern))
    if not files:
        raise SystemExit(f"no files matching {pattern} under {folder}")
    for i, f in enumerate(files, 1):
        tissue = f.name.split(".")[2]
        seen = 0
        with gzip.open(f, "rb") as fh:
            fh.readline()
            for raw in fh:
                parts = raw.split(b"\t", 3)
                if len(parts) < 3:
                    continue
                try:
                    key = (parts[0].decode(), int(parts[1]) - 1,
                           int(parts[2]) + 1)
                except ValueError:
                    continue
                if key in wanted:
                    hits[key] += 1
                    seen += 1
        print(f"  [{i}/{len(files)}] {tissue:20s} {seen:>9,} wanted junctions",
              flush=True)
    return dict(hits)


def write_csv(rows, out: Path, support_name: str,
              max_support: int | None) -> int:
    """`max_support` None leaves `confidence` empty.

    v2p copies `confidence` into every FASTA header as CONF=, so it has to
    mean support or nothing at all. Dividing by the largest value seen is
    only meaningful when the measure is graded. GTEx presence is binary -
    an event is kept exactly when all its junctions are in the file - so
    the count being divided is the event's junction arity, and a CONF of
    0.25 against 1.0 would say "weaker evidence" when it only ever meant
    "this event type has one junction, that one has four".
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["AS", "dPSI", "p_val", "Batch_support", "Lib_support",
                    "gene_id", "confidence", "event_type", support_name])
        n = 0
        for row, support in rows:
            conf = (f"{support / max_support:.4f}"
                    if max_support else "")
            w.writerow([row["event_id"], "", "", support, "", row["gene_id"],
                        conf, row["event_type"], support])
            n += 1
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True, type=Path)
    ap.add_argument("--gtex", type=Path)
    ap.add_argument("--recount3", type=Path)
    ap.add_argument("--recount3-pattern", default="*.UNIQUE.RR.gz",
                    help="UNIQUE counts uniquely-mapping reads only, the "
                         "conservative set (default); ALL includes "
                         "multi-mapping reads")
    ap.add_argument("--min-tissues", type=int, default=1,
                    help="recount3: keep a junction seen in at least this "
                         "many of the 32 tissue catalogues (default 1)")
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args(argv)

    print("reading the event catalogue ...", flush=True)
    events = list(load_events(args.events))
    usable = [(r, e, j) for r, e, j in events if j]
    wanted = {j for _, _, js in usable for j in js}
    print(f"  {len(events):,} events, {len(usable):,} with derivable "
          f"junctions, {len(wanted):,} distinct junctions needed")

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = []

    if args.gtex:
        print("scanning GTEx ...", flush=True)
        found = set(scan_gtex(args.gtex, wanted))
        kept = [(r, sum(1 for j in js if j in found))
                for r, _, js in usable if all(j in found for j in js)]
        n = write_csv(kept, args.outdir / "gtex_normal_as.csv",
                      "junctions_observed", None)
        summary.append(("GTEx v8", n, len(usable)))

    if args.recount3:
        print("scanning recount3 ...", flush=True)
        counts = scan_recount3(args.recount3, wanted, args.recount3_pattern)
        ok = {j for j, c in counts.items() if c >= args.min_tissues}
        kept = [(r, min(counts[j] for j in js))
                for r, _, js in usable if all(j in ok for j in js)]
        n = write_csv(kept, args.outdir / "recount3_normal_as.csv",
                      "min_tissues_with_junction", 32)
        summary.append(("recount3", n, len(usable)))

    print()
    for name, n, total in summary:
        print(f"{name:12s} {n:>9,} of {total:,} events "
              f"({n / total * 100:.1f}%) are observed in normal tissue")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

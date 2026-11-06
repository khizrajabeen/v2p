#!/usr/bin/env python3
"""Query splice junctions from Snaptron.

    # one region, junctions seen in at least 100 samples
    python scripts/snaptron_client.py region chr19:58351000-58353000 \
        --min-samples 100

    # every junction of a gene, written in v2p's coordinate convention
    python scripts/snaptron_client.py region BRCA1 --annotated-only \
        --out brca1_junctions.tsv

    # support for junctions listed in a file (one chr:start-end per line)
    python scripts/snaptron_client.py lookup my_junctions.txt --out support.tsv

Snaptron indexes the same junction summaries as recount, served as a
query API, so a question about a handful of genes does not require the
37 GB of count matrices.

Two conventions have to be kept apart:

    Snaptron  start = first base of the intron, end = last base
    v2p       (last base of the exon, first base of the next exon)

so a Snaptron junction becomes (start - 1, end + 1). `--convention v2p`
writes the converted pair; `--convention snaptron` leaves it alone.

Compilations are 'gtex', 'tcga', 'srav3h' and the others Snaptron hosts.
GTEx is GRCh38, verified against GENCODE v44.
"""

from __future__ import annotations

import argparse
import csv
import sys
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://snaptron.cs.jhu.edu"
TIMEOUT = 180

FIELDS = ["chromosome", "start", "end", "strand", "annotated",
          "left_motif", "right_motif", "samples_count", "coverage_sum",
          "coverage_avg"]


def query(region: str, compilation: str = "gtex",
          min_samples: int | None = None,
          min_coverage: int | None = None,
          exact: bool = False, contains: bool = False) -> list[dict]:
    """One Snaptron region query. `region` is chr:start-end or a gene name."""
    params = [("regions", region)]
    if min_samples is not None:
        params.append(("rfilter", f"samples_count>:{min_samples}"))
    if min_coverage is not None:
        params.append(("rfilter", f"coverage_sum>:{min_coverage}"))
    if exact:
        params.append(("exact", "1"))
    if contains:
        params.append(("contains", "1"))
    url = f"{BASE}/{compilation}/snaptron?" + urllib.parse.urlencode(params)

    with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:
        text = fh.read().decode("utf-8", "replace")

    lines = text.splitlines()
    if not lines:
        return []
    header = lines[0].lstrip("#").split("\t")
    out = []
    for line in lines[1:]:
        if not line.strip():
            continue
        row = dict(zip(header, line.split("\t")))
        # `samples` is one entry per sample and dwarfs everything else;
        # samples_count already carries what it is wanted for.
        row.pop("samples", None)
        out.append(row)
    return out


def to_v2p(row: dict) -> tuple[int, int]:
    """Snaptron intron bounds -> v2p's (exon end, next exon start)."""
    return int(row["start"]) - 1, int(row["end"]) + 1


def _write(rows: list[dict], out: Path | None, convention: str) -> None:
    cols = FIELDS + (["v2p_donor", "v2p_acceptor"]
                     if convention == "v2p" else [])
    fh = open(out, "w", newline="", encoding="utf-8") if out else sys.stdout
    try:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            vals = [r.get(c, "") for c in FIELDS]
            if convention == "v2p":
                a, b = to_v2p(r)
                vals += [a, b]
            w.writerow(vals)
    finally:
        if out:
            fh.close()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["region", "lookup"])
    ap.add_argument("target", help="chr:start-end, a gene name, or a file of "
                                   "regions when mode is 'lookup'")
    ap.add_argument("--compilation", default="gtex")
    ap.add_argument("--min-samples", type=int, default=None,
                    help="keep junctions seen in at least this many samples")
    ap.add_argument("--min-coverage", type=int, default=None)
    ap.add_argument("--annotated-only", action="store_true",
                    help="keep only junctions the annotation already has")
    ap.add_argument("--exact", action="store_true",
                    help="match the region's bounds exactly")
    ap.add_argument("--convention", choices=["v2p", "snaptron"],
                    default="v2p")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)

    if args.mode == "region":
        rows = query(args.target, args.compilation, args.min_samples,
                     args.min_coverage, exact=args.exact)
    else:
        regions = [l.strip() for l in Path(args.target).read_text().splitlines()
                   if l.strip() and not l.startswith("#")]
        rows = []
        for i, reg in enumerate(regions, 1):
            rows += query(reg, args.compilation, args.min_samples,
                          args.min_coverage, exact=True)
            if i % 50 == 0:
                print(f"  {i}/{len(regions)} regions", file=sys.stderr,
                      flush=True)

    if args.annotated_only:
        rows = [r for r in rows if r.get("annotated") not in ("0", "", None)]

    _write(rows, args.out, args.convention)
    print(f"{len(rows)} junctions", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Distil protein FASTAs into a normal-tissue peptide set.

    python scripts/build_normal_peptides.py \
        --fasta results/NORMAL_AS_GTEx/NORMAL_AS_GTEx.target.fasta \
        --fasta results/NORMAL_AS_recount3/NORMAL_AS_recount3.target.fasta \
        --reference ref/uniprot_human_SP.fasta \
        --out data/normal_peptides_gencode_v44.txt.gz

v2p calls a peptide novel when the reference proteome does not contain
it. That is the wrong question for neoantigen work: a peptide absent from
UniProt can still be made in healthy tissue by an alternative isoform
UniProt does not list, and such a peptide is not a tumour target. This
writes the set of peptides that occur in normal tissue and are absent
from the reference, so v2p can tell the two apart.

Peptides already in the reference are left out. v2p checks those itself,
and storing them again would multiply the file for nothing.

Isoleucine and leucine are folded together throughout, as everywhere else
in v2p: they are isobaric, so a mass spectrometer cannot separate them
and a peptide differing only by I/L is not a different peptide.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.peptides import digest                           # noqa: E402


def read_fasta(path: Path):
    hdr, buf = None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf)
                hdr, buf = line.rstrip("\n"), []
            else:
                buf.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(buf)


def peptide_keys(seq: str, mc: int, lo: int, hi: int) -> set[str]:
    return {p.replace("I", "L") for p in digest(seq, mc, lo, hi)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", action="append", required=True, type=Path,
                    help="a normal-tissue protein FASTA; repeatable")
    ap.add_argument("--reference", required=True, type=Path,
                    help="reference proteome, whose peptides are excluded")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--missed-cleavages", type=int, default=2)
    ap.add_argument("--min-len", type=int, default=6)
    ap.add_argument("--max-len", type=int, default=50)
    ap.add_argument("--skip-class", action="append",
                    default=["REFERENCE", "REFERENCE_PROTEOME"],
                    help="VT= classes to ignore; the appended reference "
                         "proteome would otherwise be digested twice")
    args = ap.parse_args(argv)

    skip = {f"VT={c}" for c in args.skip_class}

    print(f"digesting reference {args.reference.name} ...", flush=True)
    ref: set[str] = set()
    for _hdr, seq in read_fasta(args.reference):
        ref |= peptide_keys(seq, args.missed_cleavages, args.min_len,
                            args.max_len)
    print(f"  {len(ref):,} reference peptides")

    normal: set[str] = set()
    for f in args.fasta:
        n_entry = before = 0
        before = len(normal)
        for hdr, seq in read_fasta(f):
            if any(s in hdr for s in skip):
                continue
            n_entry += 1
            normal |= peptide_keys(seq, args.missed_cleavages, args.min_len,
                                   args.max_len) - ref
        print(f"  {f.name}: {n_entry:,} entries, "
              f"+{len(normal) - before:,} peptides "
              f"(union {len(normal):,})", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(args.out).endswith(".gz") else open
    with opener(args.out, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# normal-tissue tryptic peptides absent from "
                 f"{args.reference.name}\n")
        fh.write(f"# missed_cleavages={args.missed_cleavages} "
                 f"min_len={args.min_len} max_len={args.max_len} "
                 f"I/L equated\n")
        for p in sorted(normal):
            fh.write(p + "\n")
    print(f"\n{len(normal):,} peptides -> {args.out} "
          f"({args.out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

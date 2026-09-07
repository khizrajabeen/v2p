#!/usr/bin/env python3
"""Stage 1 - normalise every input file into one variant manifest.

Reference-free: runs without a genome FASTA or GTF, so it can be executed
immediately to audit the input package before the translation stage.

Usage:
  python scripts/01_parse_inputs.py \
      --res      data/HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt \
      --fusion   data/HCC1395_high_confidence_Fusion_genes_all.csv \
      --as-lr    data/HCC1395_high_confidence_AS-LR_v1.csv \
      --vcf      data/high-confidence_sSNV_sIndel_v1.sort.final.combined.sort.vcf \
      --outdir   results
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.parse.inputs import (  # noqa: E402
    parse_annovar_res, parse_as_events, parse_fusions, parse_vcf,
)
from v2p.provenance import RunLogger  # noqa: E402

MANIFEST_COLS = ["variant_id", "variant_class", "source", "genes",
                 "locus", "confidence", "payload_json"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", help="somatic SNV/InDel VCF (optional)")
    ap.add_argument("--res", help="ANNOVAR multianno RNA-editing table")
    ap.add_argument("--fusion", help="fusion gene CSV")
    ap.add_argument("--as-lr", dest="as_lr", help="alternative splicing CSV")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--logdir", default="logs")
    ap.add_argument("--no-vcf-filter", action="store_true",
                    help="keep VCF records regardless of the FILTER column")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    rl = RunLogger("01_parse_inputs", args.logdir)
    rl.add_params(no_vcf_filter=args.no_vcf_filter)

    sources = []
    if args.vcf:
        rl.add_input("vcf", args.vcf)
        sources.append(("vcf", lambda: parse_vcf(
            args.vcf, require_pass=not args.no_vcf_filter, logger=rl.log)))
    if args.res:
        rl.add_input("res", args.res)
        sources.append(("res", lambda: parse_annovar_res(args.res, logger=rl.log)))
    if args.fusion:
        rl.add_input("fusion", args.fusion)
        sources.append(("fusion", lambda: parse_fusions(args.fusion, logger=rl.log)))
    if args.as_lr:
        rl.add_input("as_lr", args.as_lr)
        sources.append(("as_lr", lambda: parse_as_events(args.as_lr, logger=rl.log)))

    if not sources:
        rl.log.error("no input files supplied")
        rl.close("error")
        return 2

    manifest = outdir / "tables" / "unified_variant_manifest.tsv"
    seen: set[str] = set()
    n_dup = 0
    with open(manifest, "w", encoding="utf-8") as out:
        out.write("\t".join(MANIFEST_COLS) + "\n")
        for name, gen in sources:
            n = 0
            for rec in gen():
                vid = rec["variant_id"]
                if vid in seen:
                    n_dup += 1
                    rl.count(f"{name}.duplicate_id")
                    continue
                seen.add(vid)
                n += 1
                rl.count(f"class.{rec['variant_class']}")
                out.write("\t".join([
                    vid, rec["variant_class"], rec["source"], rec["genes"],
                    rec["locus"], str(rec["confidence"]),
                    json.dumps(rec["payload"], separators=(",", ":")),
                ]) + "\n")
            rl.count(f"source.{name}", n)
            rl.log.info("source %s contributed %d records", name, n)

    if n_dup:
        rl.log.warning("%d duplicate variant_ids collapsed", n_dup)
    rl.add_output("manifest", manifest)
    rl.close()
    print(f"\nmanifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

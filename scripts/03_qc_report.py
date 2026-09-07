#!/usr/bin/env python3
"""Stage 3 - QC and coding-impact triage on the unified manifest.

Reference-free. Produces:
  qc/summary.md                      human-readable audit
  tables/res_recoding_sites.tsv      RNA-editing sites that change a codon,
                                     taken from the ANNOVAR AAChange fields
  tables/class_counts.tsv            record counts per variant class
  tables/expected_yield.tsv          how many FASTA entries each class should
                                     produce, so stage 2 output can be checked

Usage:
  python scripts/03_qc_report.py --manifest results/tables/unified_variant_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.parse.inputs import AS_TYPES        # noqa: E402
from v2p.provenance import RunLogger         # noqa: E402
from v2p.seqops import hgvs_p_to_change      # noqa: E402

CODING_FUNC = {"exonic", "splicing", "exonic;splicing"}
RECODING_EXFUNC = {"nonsynonymous SNV", "stopgain", "stoploss",
                   "frameshift deletion", "frameshift insertion",
                   "nonframeshift deletion", "nonframeshift insertion"}


def manifest_label(manifest: str | Path, outdir: str | Path) -> str:
    """How the manifest is named inside qc/summary.md.

    summary.md ships inside the release, so anything absolute in it makes
    two runs of the same config into different directories differ - and
    that single line was the only thing that stopped the release being
    byte-identical between runs. Naming the manifest relative to the
    output directory keeps the information and drops the part that has
    nothing to do with the data.

    Forward slashes always, so a release built on Windows and one built
    on Linux read the same.
    """
    m, o = Path(manifest), Path(outdir)
    try:
        rel = Path(os.path.relpath(m.resolve(), o.resolve())).as_posix()
    except (OSError, ValueError):
        # Different drives on Windows, or a path that cannot be resolved.
        return m.name
    # A manifest outside the output directory would be named with a
    # ../.. prefix that says more about the caller's layout than about
    # the data, so name it by file instead.
    return m.name if rel.startswith("..") else rel


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--logdir", default="logs")
    args = ap.parse_args()

    out = Path(args.outdir)
    (out / "qc").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    rl = RunLogger("03_qc_report", args.logdir)
    rl.add_input("manifest", args.manifest)

    cls = Counter()
    src = Counter()
    res_func = Counter()
    res_exfunc = Counter()
    res_recoding: list[dict] = []
    as_types = Counter()
    as_conf: list[float] = []
    fusion_rows: list[dict] = []
    genes_by_class: dict[str, set[str]] = defaultdict(set)
    chrom_counter = Counter()
    bad_alleles = 0

    with open(args.manifest, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            vclass = row["variant_class"]
            cls[vclass] += 1
            src[row["source"]] += 1
            p = json.loads(row["payload_json"])
            if row["genes"]:
                for g in row["genes"].replace(",", ";").split(";"):
                    if g.strip():
                        genes_by_class[vclass].add(g.strip())

            if vclass == "RNA_EDITING":
                res_func[p.get("func_refgene", ".")] += 1
                res_exfunc[p.get("exonic_func_refgene", ".")] += 1
                chrom_counter[p.get("chrom", "?")] += 1
                if p.get("ref") + ">" + p.get("alt") not in ("A>G", "T>C"):
                    bad_alleles += 1
                if p.get("exonic_func_refgene") in RECODING_EXFUNC:
                    for a in p.get("aachange_refgene", []):
                        parsed = hgvs_p_to_change(a.get("protein", ""))
                        res_recoding.append({
                            "locus": row["locus"],
                            "ref": p.get("ref", ""), "alt": p.get("alt", ""),
                            "gene": a.get("gene", p.get("gene_refgene", "")),
                            "transcript": a.get("transcript", ""),
                            "exon": a.get("exon", ""),
                            "cdna": a.get("cdna", ""),
                            "protein_change": a.get("protein", ""),
                            "ref_aa": parsed[0] if parsed else "",
                            "aa_pos": parsed[1] if parsed else "",
                            "alt_aa": parsed[2] if parsed else "",
                            "exonic_func": p.get("exonic_func_refgene", ""),
                            "cosmic": p.get("cosmic_coding", "."),
                            "interpro_domain": p.get("interpro_domain", "."),
                        })
            elif vclass.startswith("AS_"):
                as_types[p.get("event_type", "?")] += 1
                try:
                    as_conf.append(float(row["confidence"]))
                except ValueError:
                    pass
            elif vclass == "FUSION":
                fusion_rows.append({
                    "tag": p.get("tag", ""), "gene1": p.get("gene1", ""),
                    "gene2": p.get("gene2", ""),
                    "bp1": f"{p.get('chrom1')}:{p.get('pos1')}",
                    "bp2": f"{p.get('chrom2')}:{p.get('pos2')}",
                    "intrachromosomal":
                        str(p.get("chrom1") == p.get("chrom2")),
                    "detection": p.get("detection_type", ""),
                    "validated": p.get("validated", ""),
                    "confidence": row["confidence"],
                })

    # ---- write tables --------------------------------------------------
    def dump(path: Path, rows: list[dict]) -> None:
        if not rows:
            path.write_text("")
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    res_recoding.sort(key=lambda r: (r["gene"], str(r["aa_pos"])))
    dump(out / "tables" / "res_recoding_sites.tsv", res_recoding)
    dump(out / "tables" / "fusion_events.tsv", fusion_rows)
    dump(out / "tables" / "class_counts.tsv",
         [{"variant_class": k, "n_records": v,
           "n_distinct_genes": len(genes_by_class.get(k, ()))}
          for k, v in sorted(cls.items())])

    # expected FASTA yield, so stage 2 can be checked against it
    exp = []
    for k, v in sorted(cls.items()):
        if k == "RNA_EDITING":
            n = sum(res_exfunc[e] for e in RECODING_EXFUNC)
            note = "only codon-changing sites yield a variant protein"
        elif k == "FUSION":
            n = v * 2
            note = "both partner orientations attempted per event"
        elif k.startswith("AS_"):
            n = v
            note = "1-4 isoforms per event depending on annotation match"
        else:
            n = v
            note = "one protein per affected coding transcript"
        exp.append({"variant_class": k, "n_input": v,
                    "expected_fasta_entries_min": n, "basis": note})
    dump(out / "tables" / "expected_yield.tsv", exp)

    # ---- summary -------------------------------------------------------
    total = sum(cls.values())
    L: list[str] = []
    L.append("# HCC1395 variant package - input QC\n")
    L.append(f"Manifest: `{manifest_label(args.manifest, out)}`  \n"
             f"Total records: **{total}**\n")

    L.append("## Records per variant class\n")
    L.append("| variant class | records | distinct genes |")
    L.append("|---|---:|---:|")
    for k, v in sorted(cls.items()):
        L.append(f"| {k} | {v} | {len(genes_by_class.get(k, ()))} |")

    L.append("\n## RNA editing sites\n")
    L.append(f"Allele check: {bad_alleles} site(s) are not A>G / T>C "
             f"(expected 0 for A-to-I editing).\n")
    L.append("| genomic region | sites |")
    L.append("|---|---:|")
    for k, v in res_func.most_common():
        L.append(f"| {k} | {v} |")
    L.append("\n| exonic consequence | sites |")
    L.append("|---|---:|")
    for k, v in res_exfunc.most_common():
        L.append(f"| {k} | {v} |")
    L.append(f"\n**{len(res_recoding)}** transcript-level recoding annotations "
             f"were recovered from the ANNOVAR `AAChange.refGene` field and "
             f"written to `tables/res_recoding_sites.tsv`. These are the RNA "
             f"editing sites that alter a codon; the remaining sites are "
             f"intronic, UTR, intergenic or synonymous and contribute no "
             f"variant protein.\n")

    L.append("## Alternative splicing events\n")
    L.append("| event type | meaning | events |")
    L.append("|---|---|---:|")
    for k, v in as_types.most_common():
        L.append(f"| {k} | {AS_TYPES.get(k, '?')} | {v} |")
    if as_conf:
        L.append(f"\nConfidence: min {min(as_conf):.3f}, "
                 f"median {sorted(as_conf)[len(as_conf)//2]:.3f}, "
                 f"max {max(as_conf):.3f}\n")

    L.append("## Fusion events\n")
    intra = sum(1 for f in fusion_rows if f["intrachromosomal"] == "True")
    L.append(f"{len(fusion_rows)} events: {intra} intra-chromosomal, "
             f"{len(fusion_rows)-intra} inter-chromosomal.\n")

    (out / "qc" / "summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    for k, v in cls.items():
        rl.count(f"class.{k}", v)
    rl.count("res_recoding_annotations", len(res_recoding))
    rl.add_output("summary", out / "qc" / "summary.md")
    rl.add_output("res_recoding", out / "tables" / "res_recoding_sites.tsv")
    rl.add_output("expected_yield", out / "tables" / "expected_yield.tsv")
    rl.close()
    print((out / "qc" / "summary.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

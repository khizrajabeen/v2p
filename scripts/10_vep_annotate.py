#!/usr/bin/env python3
"""Stage 10 - annotate a VCF with Ensembl VEP over the REST API.

Two things need this, and neither could be done without it:

  1. pypgatk's `vcf-to-proteindb` reads the transcript and consequence
     from a VEP annotation field in INFO. Given a sites-only VCF it finds
     nothing, writes nothing, and exits 0 - so the cross-tool benchmark
     cannot run at all until the input is annotated.
  2. An independent second opinion on our own consequence calls
     (docs/BUILD_SPEC.md, milestone 7d). VEP is not ground truth, but a
     disagreement is worth looking at, and agreement between two
     independent implementations is worth more than either alone.

The REST API is used rather than a local VEP install because the offline
cache is roughly 25 GB. The trade is throughput: REST is rate-limited and
batched, which is fine for a few thousand variants and wrong for a whole
genome. For that, install VEP.

**This needs network.** No test imports it; the suite stays offline.

    python3 scripts/10_vep_annotate.py --vcf in.vcf --out-vcf annotated.vcf \\
        --out-tsv vep_consequences.tsv
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.provenance import RunLogger                       # noqa: E402

ENDPOINT = "https://rest.ensembl.org/vep/{species}/region"

# The CSQ subfields we emit, in this order. pypgatk is told these names
# explicitly rather than left to guess at a VEP version's default layout.
CSQ_FIELDS = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
              "Feature_type", "Feature", "BIOTYPE", "HGVSp"]

# VEP's own severity ordering, most severe first.
SEVERITY = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant",
    "TFBS_ablation", "TF_binding_site_variant", "regulatory_region_variant",
    "intergenic_variant",
]
_RANK = {c: i for i, c in enumerate(SEVERITY)}


def read_vcf(path: Path):
    """Return (header_lines, records) where a record is the split row."""
    head, rows = [], []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                head.append(line)
            elif line.strip():
                rows.append(line.split("\t"))
    return head, rows


def post(species: str, payload: list[str], retries: int = 4) -> list:
    body = json.dumps({"variants": payload}).encode()
    req = urllib.request.Request(
        ENDPOINT.format(species=species), data=body,
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"})
    delay = 2.0
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            # 429 is the documented rate limit; back off rather than
            # hammering a public service.
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < retries - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise
    return []


def most_severe(terms) -> str:
    terms = list(terms or [])
    return sorted(terms, key=lambda t: _RANK.get(t, 999))[0] if terms else ""


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--out-tsv", default="")
    ap.add_argument("--species", default="human")
    ap.add_argument("--batch", type=int, default=200,
                    help="variants per request; the REST API caps this")
    ap.add_argument("--logdir", default="logs")
    a = ap.parse_args()

    rl = RunLogger("10_vep_annotate", a.logdir)
    rl.add_params(species=a.species, batch=a.batch, endpoint=ENDPOINT)
    rl.add_input("vcf", a.vcf)

    head, rows = read_vcf(Path(a.vcf))
    rl.log.info("%d variants to annotate", len(rows))

    # VEP's region format wants the bare contig name.
    queries, keys = [], []
    for r in rows:
        chrom = r[0][3:] if r[0].startswith("chr") else r[0]
        queries.append(f"{chrom} {r[1]} . {r[3]} {r[4]} . . .")
        keys.append((r[0], r[1], r[3], r[4]))

    ann: dict[tuple, list] = {}
    for i in range(0, len(queries), a.batch):
        chunk = queries[i:i + a.batch]
        rl.log.info("batch %d-%d of %d", i + 1, i + len(chunk), len(queries))
        try:
            out = post(a.species, chunk)
        except Exception as e:                       # noqa: BLE001
            rl.log.error("batch failed: %s", e)
            rl.close(status="error")
            print(f"VEP request failed: {e}", file=sys.stderr)
            return 1
        for item in out:
            parts = (item.get("input") or "").split()
            if len(parts) < 5:
                continue
            c, pos, _id, ref, alt = parts[:5]
            k = (c if c.startswith("chr") else "chr" + c, pos, ref, alt)
            ann[k] = item.get("transcript_consequences") or []
        time.sleep(0.4)          # be polite to a public endpoint

    n_ann = sum(1 for k in keys if ann.get(k))
    rl.count("variants", len(keys))
    rl.count("annotated", n_ann)
    rl.log.info("%d of %d variants got at least one transcript consequence",
                n_ann, len(keys))

    # ---- annotated VCF, for pypgatk ------------------------------------
    out_head = [h for h in head if not h.startswith("#CHROM")]
    out_head.append(
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence '
        'annotations from Ensembl VEP REST. Format: '
        + "|".join(CSQ_FIELDS) + '">')
    out_head += [h for h in head if h.startswith("#CHROM")]

    with open(a.out_vcf, "w", newline="\n", encoding="utf-8") as fh:
        for h in out_head:
            fh.write(h + "\n")
        for r, k in zip(rows, keys):
            csq = []
            for tc in ann.get(k, []):
                csq.append("|".join([
                    tc.get("variant_allele", ""),
                    "&".join(tc.get("consequence_terms", [])),
                    tc.get("impact", ""),
                    tc.get("gene_symbol", ""),
                    tc.get("gene_id", ""),
                    "Transcript",
                    tc.get("transcript_id", ""),
                    tc.get("biotype", ""),
                    tc.get("hgvsp", ""),
                ]))
            row = list(r)
            if csq:
                info = row[7] if len(row) > 7 and row[7] != "." else ""
                row[7] = (info + ";" if info else "") + "CSQ=" + ",".join(csq)
            fh.write("\t".join(row) + "\n")
    rl.add_output("annotated_vcf", a.out_vcf)
    print(f"wrote {a.out_vcf}  ({n_ann}/{len(keys)} annotated)")

    # ---- per-variant table, for the consequence cross-check -------------
    if a.out_tsv:
        with open(a.out_tsv, "w", newline="\n", encoding="utf-8") as fh:
            fh.write("chrom\tpos\tref\talt\tvep_consequence\tvep_impact\t"
                     "vep_transcript\tvep_gene\tvep_hgvsp\n")
            for k in keys:
                tcs = ann.get(k, [])
                coding = [t for t in tcs
                          if t.get("biotype") == "protein_coding"] or tcs
                if not coding:
                    fh.write("\t".join(k) + "\t\t\t\t\t\n")
                    continue
                best = sorted(
                    coding,
                    key=lambda t: _RANK.get(
                        most_severe(t.get("consequence_terms")), 999))[0]
                fh.write("\t".join([
                    *k,
                    most_severe(best.get("consequence_terms")),
                    best.get("impact", ""),
                    best.get("transcript_id", ""),
                    best.get("gene_symbol", ""),
                    best.get("hgvsp", ""),
                ]) + "\n")
        rl.add_output("consequence_table", a.out_tsv)
        print(f"wrote {a.out_tsv}")

    rl.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

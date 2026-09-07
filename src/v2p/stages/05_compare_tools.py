#!/usr/bin/env python3
"""Stage 5 - compare protein FASTAs produced by different tools.

Different tools name entries differently (VEP uses Ensembl protein ids,
AGFusion uses its own scheme, ours uses UniProt-style ids), so keyed
comparison alone under-reports agreement. Two complementary views:

  1. SEQUENCE-SET overlap - how many complete protein sequences are
     identical between tools. Strict; sensitive to transcript choice.
  2. K-MER SPACE overlap - the set of distinct k-mers (default 9, the
     MHC-I / tryptic-peptide scale). This is what actually determines
     whether a downstream search or binding prediction can find a
     peptide. Two tools can pick different transcripts and still cover
     the same k-mer space, which is the outcome that matters.

Optionally restricts to variant-proximal k-mers only, since agreement on
shared reference regions is uninformative.

Usage:
  python src/v2p/stages/05_compare_tools.py \
      --fasta ours=results/fasta/HCC1395_variant_proteins.uniprot.fasta \
      --fasta vep=vep_out/mutated.fa \
      --fasta pgatk=pgatk_out/variants.fa \
      --k 9 --outdir results
"""

from __future__ import annotations

import argparse
import csv
import itertools
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from v2p.provenance import RunLogger      # noqa: E402

_GENE_PATTERNS = [
    re.compile(r"\bGN=([A-Za-z0-9._:-]+)"),        # UniProt / ours
    re.compile(r"\\GName=([A-Za-z0-9._:-]+)"),     # PEFF
    re.compile(r"\bgene_symbol:([A-Za-z0-9._-]+)"),  # Ensembl FASTA dumps
    re.compile(r"\bsymbol=([A-Za-z0-9._-]+)"),
]


def read_fasta(path: Path):
    name, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def guess_gene(header: str) -> str:
    for pat in _GENE_PATTERNS:
        m = pat.search(header)
        if m:
            return m.group(1)
    # AGFusion / pVACfuse style: GENE1-GENE2_...  or  MT.GENE.TX...
    parts = re.split(r"[|\s.]", header)
    for p in parts:
        if re.fullmatch(r"[A-Z0-9]{2,}(--?[A-Z0-9]{2,})?", p) and not p.isdigit():
            return p
    return "NA"


def kmers(seq: str, k: int) -> set[str]:
    s = seq.replace("*", "")
    return {s[i:i + k] for i in range(len(s) - k + 1)} if len(s) >= k else set()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", action="append", required=True,
                    metavar="NAME=PATH",
                    help="repeatable; label the tool, e.g. ours=out.fasta")
    ap.add_argument("--k", type=int, default=9)
    ap.add_argument("--min-len", type=int, default=8)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--logdir", default="logs")
    args = ap.parse_args()

    out = Path(args.outdir)
    (out / "qc").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    rl = RunLogger("05_compare_tools", args.logdir)
    rl.add_params(k=args.k, min_len=args.min_len)

    sets: dict[str, dict] = {}
    for spec in args.fasta:
        if "=" not in spec:
            rl.log.error("--fasta needs NAME=PATH, got %r", spec)
            rl.close("error")
            return 2
        name, path = spec.split("=", 1)
        rl.add_input(name, path)
        seqs: set[str] = set()
        km: set[str] = set()
        genes: Counter = Counter()
        by_gene: dict[str, set[str]] = {}
        n = 0
        for hdr, seq in read_fasta(Path(path)):
            if len(seq) < args.min_len:
                continue
            n += 1
            seqs.add(seq)
            km |= kmers(seq, args.k)
            g = guess_gene(hdr)
            genes[g] += 1
            by_gene.setdefault(g, set()).add(seq)
        sets[name] = {"seqs": seqs, "kmers": km, "genes": genes,
                      "by_gene": by_gene, "n_entries": n}
        rl.log.info("%s: %d entries, %d distinct sequences, %d distinct %d-mers, "
                    "%d genes", name, n, len(seqs), len(km), args.k, len(genes))
        rl.count(f"{name}.entries", n)

    names = list(sets)
    L = ["# Tool comparison\n",
         f"k-mer size **{args.k}**, minimum protein length {args.min_len}.\n"]
    L.append("| tool | entries | distinct sequences | distinct k-mers | genes |")
    L.append("|---|---:|---:|---:|---:|")
    for nme in names:
        s = sets[nme]
        L.append(f"| {nme} | {s['n_entries']} | {len(s['seqs'])} | "
                 f"{len(s['kmers'])} | {len(s['genes'])} |")

    rows = []
    L.append("\n## Pairwise agreement\n")
    L.append("| A | B | shared sequences | Jaccard (seq) | shared k-mers | "
             "Jaccard (k-mer) | A-only genes | B-only genes |")
    L.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for a, b in itertools.combinations(names, 2):
        sa, sb = sets[a], sets[b]
        si = sa["seqs"] & sb["seqs"]
        su = sa["seqs"] | sb["seqs"]
        ki = sa["kmers"] & sb["kmers"]
        ku = sa["kmers"] | sb["kmers"]
        ga = set(sa["genes"]) - set(sb["genes"])
        gb = set(sb["genes"]) - set(sa["genes"])
        js = len(si) / len(su) if su else 0.0
        jk = len(ki) / len(ku) if ku else 0.0
        L.append(f"| {a} | {b} | {len(si)} | {js:.3f} | {len(ki)} | "
                 f"{jk:.3f} | {len(ga)} | {len(gb)} |")
        rows.append({"tool_a": a, "tool_b": b, "shared_sequences": len(si),
                     "jaccard_sequence": round(js, 4),
                     "shared_kmers": len(ki), "jaccard_kmer": round(jk, 4),
                     "genes_only_in_a": ";".join(sorted(ga)[:50]),
                     "genes_only_in_b": ";".join(sorted(gb)[:50])})
        rl.count(f"pair.{a}_vs_{b}.shared_sequences", len(si))

        # genes present in both but with no identical sequence: the
        # informative disagreements, where both tools translated the same
        # gene and got different answers
        both = set(sa["genes"]) & set(sb["genes"])
        disagree = sorted(g for g in both
                          if not (sa["by_gene"][g] & sb["by_gene"][g]))
        if disagree:
            L.append(f"\n<details><summary>{len(disagree)} genes translated by "
                     f"both {a} and {b} with no identical sequence</summary>\n")
            L.append("`" + " ".join(disagree[:200]) + "`\n")
            L.append("</details>\n")
            rl.count(f"pair.{a}_vs_{b}.genes_disagree", len(disagree))

    if rows:
        with open(out / "tables" / "tool_comparison.tsv", "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    L.append("\n## How to read this\n")
    L.append("- High k-mer Jaccard with low sequence Jaccard means the tools "
             "agree on the peptides but chose different transcripts. That is "
             "usually acceptable — decide the transcript policy and move on.\n")
    L.append("- Low k-mer Jaccard means a real disagreement in translation. "
             "Inspect the listed genes individually before shipping.\n")
    L.append("- Genes present in only one tool are usually a scope difference "
             "(e.g. VEP does not do fusions) rather than an error.\n")

    (out / "qc" / "tool_comparison.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")
    rl.add_output("report", out / "qc" / "tool_comparison.md")
    rl.close()
    print((out / "qc" / "tool_comparison.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

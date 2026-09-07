#!/usr/bin/env python3
"""Stage 4 - validate against the UniProt reference proteome.

Two independent checks:

  A. ANNOTATION CHECK (no reference genome needed). For every variant that
     carries a reported protein change, confirm the *reference* residue at
     that position really is what the annotation claims, according to
     UniProt. Catches transcript/isoform mismatches and wrong-build calls.

  B. TRANSLATION CHECK (needs stage 2 output with --include-reference).
     Compare our translated reference proteins against UniProt residue by
     residue. If our reference GRIA2 does not equal UniProt's GRIA2, our
     variant GRIA2 is not trustworthy either. This is the strongest QC
     gate available and it needs no second tool.

Usage:
  # check A only, runs today
  python scripts/04_validate_uniprot.py --uniprot ref/uniprot_human_SP.fasta \
      --recoding results/tables/res_recoding_sites.tsv

  # A + B, after stage 2
  python scripts/04_validate_uniprot.py --uniprot ref/uniprot_human_SP.fasta \
      --recoding results/tables/res_recoding_sites.tsv \
      --protein-fasta results/fasta/HCC1395_variant_proteins.uniprot.fasta
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from v2p.provenance import RunLogger        # noqa: E402
from v2p.validate import (                  # noqa: E402
    UniProtDB, best_uniprot_match, check_reference_residue,
    compare_sequences, selenoprotein_genes,
)

# Statuses that mean "the same protein", as opposed to a translation error.
# UniProt canonical and GENCODE canonical are chosen independently and
# frequently disagree on the start codon, so an N-terminal offset is an
# annotation-source difference, not a defect in the conversion.
AGREEING = {"identical", "identical_after_met_trim",
            "isoform_extension", "offset_isoform",
            "selenoprotein_expected_truncation"}

# fail the run if the reference-translation identity rate drops below this
DEFAULT_MIN_IDENTITY_RATE = 0.90


def read_gencode_translations(path: Path) -> dict[str, str]:
    """Map ENST id -> protein, from a GENCODE pc_translations FASTA.

    Header layout: >ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|name|symbol|length
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, str] = {}
    name, buf = None, []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).rstrip("*")
                buf = []
                parts = line[1:].split("|")
                name = parts[1] if len(parts) > 1 else None
            elif name:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).rstrip("*")
    return out


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
                name, chunks = line, []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parse_kv(header: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for tok in header.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--recoding", help="res_recoding_sites.tsv from stage 3")
    ap.add_argument("--protein-fasta", help="stage 2 FASTA (with REFERENCE entries)")
    ap.add_argument("--gencode-translations",
                    help="gencode.vNN.pc_translations.fa.gz - GENCODE's own "
                         "translations of the same transcript ids. This is "
                         "the authoritative correctness check: it isolates "
                         "our translation from the choice of canonical "
                         "isoform, which UniProt comparison cannot do.")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--logdir", default="logs")
    ap.add_argument("--min-identity-rate", type=float,
                    default=DEFAULT_MIN_IDENTITY_RATE)
    args = ap.parse_args()

    out = Path(args.outdir)
    (out / "qc").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    rl = RunLogger("04_validate_uniprot", args.logdir)
    rl.add_params(min_identity_rate=args.min_identity_rate)
    rl.add_input("uniprot", args.uniprot)

    db = UniProtDB.load(args.uniprot, logger=rl.log)
    seleno = selenoprotein_genes(db)
    rl.log.info("selenoprotein genes in the reference: %d (now translated "
                "with U via the GENCODE Selenocysteine annotation, so they "
                "are compared normally rather than excluded)", len(seleno))
    if db.malformed_headers:
        rl.log.warning("%d malformed UniProt headers", len(db.malformed_headers))

    L: list[str] = ["# UniProt validation report\n"]
    L.append(f"Reference: `{Path(args.uniprot).name}` — "
             f"{len(db.entries)} entries, {len(db.by_gene)} genes.\n")
    status = "ok"

    # ---------------- check A -------------------------------------------
    if args.recoding:
        rl.add_input("recoding", args.recoding)
        rows = list(csv.DictReader(open(args.recoding, encoding="utf-8"),
                                   delimiter="\t"))
        per_site: dict[tuple, list] = defaultdict(list)
        detail = []
        st = Counter()
        for r in rows:
            if not r.get("aa_pos"):
                continue
            rc = check_reference_residue(db, r["gene"], int(r["aa_pos"]),
                                         r["ref_aa"])
            st[rc.status] += 1
            per_site[(r["locus"], r["gene"])].append(rc.status)
            detail.append({
                "locus": r["locus"], "gene": r["gene"],
                "transcript": r["transcript"],
                "protein_change": r["protein_change"],
                "annotated_ref_aa": r["ref_aa"],
                "uniprot_aa": rc.observed_aa,
                "uniprot_accession": rc.accession,
                "uniprot_length": rc.uniprot_length,
                "status": rc.status,
            })
        with open(out / "tables" / "uniprot_residue_check.tsv", "w",
                  newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(detail[0]), delimiter="\t")
            w.writeheader()
            w.writerows(detail)

        anchored = sum(1 for v in per_site.values() if "match" in v)
        absent = sum(1 for v in per_site.values() if set(v) == {"no_uniprot_entry"})
        mism = len(per_site) - anchored - absent
        for k, v in st.items():
            rl.count(f"residue_check.{k}", v)
        rl.count("sites.anchored", anchored)
        rl.count("sites.all_mismatch", mism)

        L.append("## A. Annotation check — reference residue vs UniProt\n")
        L.append(f"{len(detail)} transcript-level annotations across "
                 f"{len(per_site)} unique sites.\n")
        L.append("| outcome | annotations |")
        L.append("|---|---:|")
        for k, v in st.most_common():
            L.append(f"| {k} | {v} |")
        L.append(f"\nPer site: **{anchored}** anchored (at least one "
                 f"transcript's reference residue agrees with UniProt), "
                 f"{absent} gene absent from SwissProt, {mism} with no "
                 f"agreeing transcript.\n")
        L.append("A mismatch is normally an isoform-numbering difference "
                 "between the annotation's RefSeq transcript and the UniProt "
                 "canonical entry, not a wrong call — but each one is listed "
                 "in `tables/uniprot_residue_check.tsv` so it can be checked "
                 "rather than assumed.\n")
        if anchored == 0 and len(per_site):
            status = "error"
            rl.log.error("no site anchored to UniProt — check the reference build")

    # ---------------- check B ---------------------------------------
    if args.protein_fasta:
        rl.add_input("protein_fasta", args.protein_fasta)
        refs: dict[str, tuple[str, str]] = {}   # tx -> (gene, sequence)
        for hdr, seq in read_fasta(Path(args.protein_fasta)):
            kv = parse_kv(hdr)
            vt = kv.get("VT") or kv.get("VariantType", "")
            if "REFERENCE" not in (vt + hdr):
                continue
            gene = kv.get("GN") or kv.get("GName", "")
            tx = kv.get("TX") or kv.get("TranscriptId", "")
            if gene and seq:
                refs[tx or f"nolabel:{gene}:{len(refs)}"] = (gene, seq)
        rl.log.info("%d reference proteins recovered from the FASTA", len(refs))

        # ---- B1: against GENCODE's own translations (authoritative) ----
        gencode_rate = None
        if args.gencode_translations:
            rl.add_input("gencode_translations", args.gencode_translations)
            gc = read_gencode_translations(Path(args.gencode_translations))
            rl.log.info("GENCODE translations: %d transcripts", len(gc))
            g_st = Counter()
            g_bad = []
            for tx, (gene, seq) in refs.items():
                theirs = gc.get(tx)
                if theirs is None:
                    g_st["transcript_not_in_gencode"] += 1
                    continue
                c = compare_sequences(seq, theirs)
                g_st[c["status"]] += 1
                if c["status"] != "identical":
                    g_bad.append({"transcript": tx, "gene": gene,
                                  "our_length": c["len_ours"],
                                  "gencode_length": c["len_theirs"],
                                  "identity": c["identity"],
                                  "status": c["status"]})
            n_cmp = sum(v for k, v in g_st.items()
                        if k != "transcript_not_in_gencode")
            gencode_rate = g_st["identical"] / n_cmp if n_cmp else None
            for k, v in g_st.items():
                rl.count(f"gencode_check.{k}", v)
            if g_bad:
                with open(out / "tables" / "gencode_translation_check.tsv",
                          "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(g_bad[0]),
                                       delimiter="\t")
                    w.writeheader()
                    w.writerows(g_bad)

            L.append("\n## B1. Translation check vs GENCODE (authoritative)\n")
            L.append(f"{n_cmp} transcripts compared against GENCODE's own "
                     f"translation of the same transcript ids. "
                     f"Exact-agreement: **{gencode_rate:.1%}**.\n"
                     if gencode_rate is not None else
                     "No transcripts could be compared.\n")
            L.append("| outcome | transcripts |")
            L.append("|---|---:|")
            for k, v in g_st.most_common():
                L.append(f"| {k} | {v} |")
            L.append("\nThis is the check that decides whether the conversion "
                     "is correct. It uses the *same* transcript ids we "
                     "translated, so it isolates exon assembly, CDS offset, "
                     "strand handling and the codon table from the separate "
                     "question of which isoform is canonical.\n")
            if gencode_rate is not None and gencode_rate < args.min_identity_rate:
                status = "error"
                rl.log.error("GENCODE agreement %.1f%% is below the %.0f%% "
                             "threshold", gencode_rate * 100,
                             args.min_identity_rate * 100)

        # ---- B2: against UniProt (informational) ------------------------
        cmp_rows = []
        st2 = Counter()
        for tx, (gene, seq) in refs.items():
            up, c = best_uniprot_match(db, gene, seq)
            if up is None:
                st2["no_uniprot_entry"] += 1
                continue
            st2[c["status"]] += 1
            cmp_rows.append({
                "gene": gene, "transcript": tx,
                "uniprot_accession": up.accession,
                "our_length": c["len_ours"], "uniprot_length": c["len_theirs"],
                "offset": c.get("offset", 0), "identity": c["identity"],
                "first_diff": c["first_diff"], "status": c["status"],
            })
        if cmp_rows:
            with open(out / "tables" / "uniprot_translation_check.tsv", "w",
                      newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(cmp_rows[0]),
                                   delimiter="\t")
                w.writeheader()
                w.writerows(cmp_rows)

        n_ref = sum(st2.values())
        agree = sum(v for k, v in st2.items() if k in AGREEING)
        rate = agree / n_ref if n_ref else 0.0
        for k, v in st2.items():
            rl.count(f"uniprot_translation_check.{k}", v)
        L.append("\n## B2. Comparison vs UniProt canonical (informational)\n")
        L.append(f"{n_ref} reference proteins. Same-protein rate "
                 f"(identical, Met-trim, isoform extension or fixed offset): "
                 f"**{rate:.1%}**.\n")
        L.append("| outcome | proteins |")
        L.append("|---|---:|")
        for k, v in st2.most_common():
            L.append(f"| {k} | {v} |")
        L.append("\nUniProt canonical and GENCODE canonical are curated "
                 "independently and often pick different start codons, so a "
                 "residue offset here is an annotation-source difference "
                 "rather than a conversion error. Only `divergent` entries "
                 "warrant inspection; they are listed in "
                 "`tables/uniprot_translation_check.tsv`.\n")
        if gencode_rate is None and n_ref and rate < args.min_identity_rate:
            status = "error"
            rl.log.error("UniProt same-protein rate %.1f%% is below the %.0f%% "
                         "threshold and no GENCODE file was supplied",
                         rate * 100, args.min_identity_rate * 100)

    (out / "qc" / "uniprot_validation.md").write_text("\n".join(L) + "\n",
                                                      encoding="utf-8")
    rl.add_output("report", out / "qc" / "uniprot_validation.md")
    rl.close(status)
    print((out / "qc" / "uniprot_validation.md").read_text())
    return 1 if status == "error" else 0


if __name__ == "__main__":
    raise SystemExit(main())

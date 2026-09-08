#!/usr/bin/env python3
"""Stage 2 - translate the unified manifest into a combined protein FASTA.

Requires a reference genome FASTA (hg38/GRCh38, same build as the input
calls) and a matching GENCODE/Ensembl GTF.

Usage:
  python src/v2p/stages/02_build_protein_fasta.py \
      --manifest results/tables/unified_variant_manifest.tsv \
      --genome   ref/GRCh38.primary_assembly.genome.fa \
      --gtf      ref/gencode.v44.annotation.gtf.gz \
      --outdir   results \
      --header-style peff
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from v2p.annotation import Annotation, Genome            # noqa: E402
from v2p.build.fusion import build_fusion_proteins       # noqa: E402
from v2p.build.smallvar import (                         # noqa: E402
    ProteinRecord, build_small_variant_proteins,
)
from v2p.build.splicing import build_splicing_proteins   # noqa: E402
from v2p.build.noncanonical import (
    build_noncanonical_proteins, build_utr_orfs, drop_known_proteins,
)
from v2p.species import load_species
from v2p.fasta import (                                  # noqa: E402
    HEADER_STYLES, write_fasta, write_record_table,
)
from v2p.provenance import RunLogger                     # noqa: E402
from v2p.validate import UniProtDB, check_reference_residue  # noqa: E402

SMALL_CLASSES = {"SNV", "MNV", "INDEL", "RNA_EDITING"}

# Consequences that leave the protein byte-identical to the reference.
# They are real variants, but they contribute no sequence to a search
# database, and they are three different biological situations that
# should not be reported as one.
NO_PROTEIN_CHANGE = {"synonymous", "5_prime_UTR", "3_prime_UTR"}


def load_manifest(path: Path):
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["payload"] = json.loads(row["payload_json"])
            yield row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--genome", required=True, help="reference FASTA (bgzip or plain)")
    ap.add_argument("--gtf", required=True, help="GENCODE/Ensembl GTF")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--logdir", default="logs")
    ap.add_argument("--header-style", default="descriptive",
                    choices=sorted(HEADER_STYLES))
    ap.add_argument("--transcript-mode", default="representative",
                    choices=["representative", "all"])
    ap.add_argument("--species", default="human",
                    help="species name (human, mouse, ...) or a path to a "
                         "config/species/*.yaml. Sets the entry-name suffix "
                         "and the OS=/OX= fields in headers.")
    ap.add_argument("--include-noncanonical", action="store_true",
                    help="also three-frame translate non-coding "
                         "transcripts (lncRNA, pseudogene). OFF by "
                         "default: it multiplies database size and an "
                         "inflated search space costs sensitivity.")
    ap.add_argument("--nc-min-aa", type=int, default=30,
                    help="minimum ORF length in residues (default 30)")
    ap.add_argument("--nc-any-start", action="store_true",
                    help="keep ORFs that do not begin at ATG, recording "
                         "the frame")
    ap.add_argument("--classes", default="",
                    help="comma-separated variant_class filter (default: all)")
    ap.add_argument("--min-peptide", type=int, default=8,
                    help="drop proteins shorter than this many residues")
    ap.add_argument("--keep-synonymous", action="store_true",
                    help="keep variants whose protein is identical to the "
                         "reference. Off by default: a synonymous change "
                         "yields a byte-identical sequence, so the entry adds "
                         "no peptide and only inflates the search space.")
    ap.add_argument("--include-reference", action="store_true",
                    help="also emit the unmodified reference protein of each "
                         "affected transcript (recommended for MS search)")
    ap.add_argument("--no-construct-isoforms", action="store_true",
                    help="splicing: annotated-transcript matching only")
    ap.add_argument("--uniprot", default="",
                    help="UniProt SwissProt FASTA; enables per-record "
                         "accession/description lookup and residue QC")
    ap.add_argument("--append-full-proteome", action="store_true",
                    help="append every UniProt entry verbatim, so the file is "
                         "a drop-in replacement for the reference database")
    ap.add_argument("--emit-disposition", default="",
                    help="write one row per input variant recording whether "
                         "it produced a protein and, if not, why. This is the "
                         "only way to state a recovery rate honestly: a "
                         "variant that yields nothing leaves no trace in the "
                         "FASTA.")
    ap.add_argument("--emit-isoform-gtf", default="",
                    help="write the splicing isoform structures to a GTF, for "
                         "independent translation with `gffread -y`")
    ap.add_argument("--limit", type=int, default=0, help="debug: cap records")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    rl = RunLogger("02_build_protein_fasta", args.logdir)
    species = load_species(args.species)
    rl.add_params(species=species.common_name,
                  taxon_id=species.taxon_id,
                  header_style=args.header_style,
                  transcript_mode=args.transcript_mode,
                  min_peptide=args.min_peptide,
                  include_reference=args.include_reference,
                  keep_synonymous=args.keep_synonymous,
                  construct_isoforms=not args.no_construct_isoforms,
                  classes=args.classes or "ALL")
    rl.add_input("manifest", args.manifest)
    rl.add_input("genome", args.genome, checksum=False)
    rl.add_input("gtf", args.gtf)

    rl.log.info("loading annotation ...")
    # Non-coding transcripts are dropped at load time unless the
    # three-frame path needs them; keeping them otherwise would slow
    # every ordinary run for nothing.
    ann = Annotation.from_gtf(args.gtf,
                              coding_only=not args.include_noncanonical,
                              logger=rl.log)
    ann.build_position_index(logger=rl.log)
    rl.log.info("opening genome ...")
    genome = Genome(args.genome)

    updb = None
    if args.uniprot:
        rl.add_input("uniprot", args.uniprot)
        updb = UniProtDB.load(args.uniprot, logger=rl.log)

    wanted = {c.strip() for c in args.classes.split(",") if c.strip()}
    records: list[ProteinRecord] = []
    disposition: list[dict] = []
    ref_emitted: set[str] = set()
    n_rows = 0

    for row in load_manifest(Path(args.manifest)):
        vclass = row["variant_class"]
        if wanted and vclass not in wanted:
            continue
        n_rows += 1
        if args.limit and n_rows > args.limit:
            break
        p = row["payload"]
        rl.count(f"input.{vclass}")
        try:
            if vclass in SMALL_CLASSES:
                genes = p.get("gene_refgene") or row["genes"]
                recs = build_small_variant_proteins(
                    p["chrom"], p["pos"], p["ref"], p["alt"], ann, genome,
                    variant_class=vclass, genes=genes, source=row["source"],
                    confidence=row["confidence"],
                    transcript_mode=args.transcript_mode, logger=rl.log)
            elif vclass == "FUSION":
                recs = build_fusion_proteins(
                    p["gene1"], p["chrom1"], p["pos1"],
                    p["gene2"], p["chrom2"], p["pos2"], ann, genome,
                    source=row["source"], confidence=row["confidence"],
                    detection_type=p.get("detection_type", ""),
                    validated=p.get("validated", ""), logger=rl.log)
            elif vclass.startswith("AS_"):
                recs = build_splicing_proteins(
                    p["event_id"], p["event_type"], p["gene_id"], p["chrom"],
                    p["strand"], p["coord_groups"], ann, genome,
                    source=row["source"], confidence=row["confidence"],
                    dpsi=p.get("dPSI"), p_val=p.get("p_val"),
                    allow_construction=not args.no_construct_isoforms,
                    logger=rl.log)
            else:
                rl.count(f"skipped.unknown_class.{vclass}")
                continue
        except Exception as exc:                      # keep the run going
            rl.count(f"error.{vclass}")
            rl.log.error("%s %s failed: %s: %s", vclass, row["variant_id"],
                         type(exc).__name__, exc)
            disposition.append({"variant_id": row["variant_id"],
                                "variant_class": vclass, "n_proteins": 0,
                                "outcome": "error",
                                "detail": type(exc).__name__})
            continue

        n_before = len(recs)
        consequences = sorted({r.consequence for r in recs})

        n_kept_here = 0
        n_syn_here = 0
        n_short_here = 0
        syn_kinds: set[str] = set()
        for r in recs:
            if len(r.sequence) < args.min_peptide:
                rl.count("dropped.too_short")
                rl.count(f"dropped.too_short.{vclass}")
                n_short_here += 1
                continue
            if (r.consequence in NO_PROTEIN_CHANGE
                    and not args.keep_synonymous):
                rl.count(f"dropped.{r.consequence}")
                rl.count(f"dropped.{r.consequence}.{vclass}")
                n_syn_here += 1
                syn_kinds.add(r.consequence)
                continue
            n_kept_here += 1
            if updb is not None:
                up = updb.longest_for_gene(r.gene.split("--")[0])
                if up is not None:
                    r.extra.setdefault("accession", up.accession)
                    r.extra.setdefault("description", up.description)
                if r.variant_pos_aa and r.extra.get("ref_protein"):
                    ref_p = r.extra["ref_protein"]
                    if r.variant_pos_aa <= len(ref_p):
                        rc = check_reference_residue(
                            updb, r.gene, r.variant_pos_aa,
                            ref_p[r.variant_pos_aa - 1])
                        r.extra["uniprot_check"] = rc.status
                        rl.count(f"uniprot_check.{rc.status}")
            rl.count(f"built.{r.variant_class}.{r.consequence}")
            records.append(r)
            if args.include_reference:
                ref_seq = r.extra.get("ref_protein")
                key = f"{r.transcript}"
                if ref_seq and key not in ref_emitted:
                    ref_emitted.add(key)
                    records.append(ProteinRecord(
                        seq_id=f"REF_{r.gene}_{r.transcript}",
                        sequence=ref_seq, variant_class="REFERENCE",
                        consequence="reference", gene=r.gene,
                        transcript=r.transcript, locus=r.locus,
                        source="reference", notes=["reference_protein"]))

        # Why did this variant yield nothing? Distinguishing "not in a
        # coding exon" from "no transcript here at all" is the difference
        # between an expected result and a lookup failure.
        if n_kept_here:
            outcome, detail = "protein_built", ";".join(consequences)
        elif n_syn_here:
            if syn_kinds <= {"5_prime_UTR", "3_prime_UTR"}:
                outcome = "utr_variant"
                detail = ";".join(sorted(syn_kinds))
            elif "synonymous" in syn_kinds and len(syn_kinds) == 1:
                outcome = "synonymous_only"
                detail = "protein identical to reference"
            else:
                outcome = "synonymous_or_utr"
                detail = ";".join(sorted(syn_kinds))
        elif n_short_here:
            outcome, detail = "below_min_length", f"<{args.min_peptide} aa"
        elif n_before:
            outcome, detail = "filtered", ";".join(consequences)
        elif vclass in SMALL_CLASSES:
            near = ann.transcripts_at(p["chrom"], p["pos"])
            if not near:
                outcome = "no_coding_transcript_at_locus"
                detail = "intergenic, or the gene has no coding transcript"
            else:
                outcome = "not_in_coding_exon"
                detail = (f"{len(near)} coding transcript(s) span this "
                          f"position but the variant is intronic, UTR, or "
                          f"crosses a splice junction")
        else:
            outcome, detail = "no_protein", "see the log for warnings"

        disposition.append({"variant_id": row["variant_id"],
                            "variant_class": vclass,
                            "n_proteins": n_kept_here,
                            "outcome": outcome, "detail": detail})
        rl.count(f"outcome.{vclass}.{outcome}")

        if n_rows % 1000 == 0:
            rl.log.info("processed %d manifest rows, %d sequences so far",
                        n_rows, len(records))

    if args.emit_disposition:
        dp = Path(args.emit_disposition)
        dp.parent.mkdir(parents=True, exist_ok=True)
        with open(dp, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["variant_id", "variant_class",
                                               "n_proteins", "outcome",
                                               "detail"],
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(disposition)
        rl.log.info("wrote disposition for %d input variants to %s",
                    len(disposition), dp)
        rl.add_output("disposition", dp)

    if args.emit_isoform_gtf:
        gpath = Path(args.emit_isoform_gtf)
        gpath.parent.mkdir(parents=True, exist_ok=True)
        n_gtf = 0
        with open(gpath, "w", encoding="utf-8") as gh:
            for r in records:
                ex = r.extra.get("gtf_exons")
                if not ex:
                    continue
                n_gtf += 1
                tid = r.seq_id
                attr = (f'gene_id "{r.gene}"; transcript_id "{tid}"; '
                        f'gene_name "{r.gene}"; variant_type "{r.variant_class}";')
                for s_, e_ in sorted(ex):
                    gh.write(f"{r.extra['gtf_chrom']}\tv2p\texon\t{s_}\t{e_}"
                             f"\t.\t{r.extra['gtf_strand']}\t.\t{attr}\n")
                for s_, e_ in sorted(r.extra.get("gtf_cds") or []):
                    gh.write(f"{r.extra['gtf_chrom']}\tv2p\tCDS\t{s_}\t{e_}"
                             f"\t.\t{r.extra['gtf_strand']}\t0\t{attr}\n")
        rl.log.info("wrote %d isoform models to %s", n_gtf, gpath)
        rl.add_output("isoform_gtf", gpath)

    if args.append_full_proteome and updb is not None:
        for acc, e in updb.entries.items():
            records.append(ProteinRecord(
                seq_id=acc, sequence=e.sequence, variant_class="REFERENCE",
                consequence="reference", gene=e.gn.split()[0] if e.gn else "",
                transcript="", source="uniprot",
                notes=[], extra={"uniprot_header": e.header(),
                                 "accession": acc}))
        rl.log.info("appended %d UniProt entries verbatim", len(updb.entries))

    tag = f"{args.header_style}.{args.transcript_mode}"
    fasta_path = outdir / "fasta" / f"HCC1395_variant_proteins.{tag}.fasta"
    table_path = outdir / "tables" / f"protein_records.{tag}.tsv"
    if args.include_noncanonical:
        rl.log.info('three-frame translating non-coding transcripts ...')
        nc = build_noncanonical_proteins(
            ann.tx.values(), genome, min_aa=args.nc_min_aa,
            require_atg=not args.nc_any_start)
        rl.log.info('non-canonical ORFs: %d', len(nc))
        rl.count('noncanonical.orfs', len(nc))
        for r in nc:
            rl.count(f'noncanonical.{r.variant_class}')
        utr = build_utr_orfs(ann.tx.values(), genome,
                             min_aa=args.nc_min_aa,
                             require_atg=not args.nc_any_start)
        rl.log.info('UTR ORFs: %d', len(utr))
        rl.count('noncanonical.NC_UTR', len(utr))
        nc_all = nc + utr
        if updb is not None:
            # A pseudogene ORF often reproduces its parent gene's protein
            # exactly. Those entries cannot yield a peptide the reference
            # does not already explain, so they inflate the search space
            # for nothing.
            nc_all, n_dup = drop_known_proteins(
                nc_all, (e.sequence for e in updb.entries.values()),
                logger=rl.log)
            rl.count("noncanonical.dropped_identical_to_reference", n_dup)
        records.extend(nc_all)

    counts = write_fasta(records, fasta_path, style=args.header_style,
                         logger=rl.log, species=species)
    write_record_table(records, table_path)
    for k, v in counts.items():
        rl.count(f"fasta.{k}", v)
    rl.add_output("fasta", fasta_path)
    rl.add_output("record_table", table_path)
    rl.close()
    print(f"\nFASTA: {fasta_path}\ntable: {table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

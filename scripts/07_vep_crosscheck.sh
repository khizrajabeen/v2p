#!/usr/bin/env bash
# Stage 7 (optional) - orthogonal cross-check with Ensembl VEP.
#
#   bash scripts/07_vep_crosscheck.sh <vcf> <res_multianno> <genome> <vep_cache_dir>
#
# Not required for correctness: agreement with GENCODE's own translations
# is already 100%, which establishes that the conversion is right. This
# adds an independently implemented second opinion on *consequence
# classification*, which is what a reviewer or a methods section may want
# named. It costs about an hour of single-threaded CPU.

set -euo pipefail

VCF="${1:?usage: $0 <vcf> <res_multianno> <genome.fa> <vep_cache_dir>}"
RES="${2:?}"
GENOME="${3:?}"
CACHE="${4:-$HOME/.vep}"
OUT="vep_out"
mkdir -p "$OUT"

command -v vep >/dev/null || {
  echo "vep not found. Install with:"
  echo "  conda install -c bioconda ensembl-vep=112"
  echo "  vep_install -a cfp -s homo_sapiens -y GRCh38 -c $CACHE -g ProteinSeqs,Wildtype,Frameshift"
  exit 1
}

echo "== 1/3 VEP on somatic SNV/InDel =="
# --pick gives one transcript per variant, matching our default
# --transcript-mode representative. Drop it and use --transcript-mode all
# on our side if you want a like-for-like all-transcript comparison.
# ProteinSeqs disables --fork, so this is single-threaded by design.
vep --input_file "$VCF" --output_file "$OUT/snv_indel.txt" --force_overwrite \
    --species homo_sapiens --assembly GRCh38 --cache --offline \
    --dir_cache "$CACHE" --fasta "$GENOME" \
    --symbol --protein --uniprot --hgvs --coding_only --pick \
    --plugin ProteinSeqs,"$OUT/snv_reference.fa","$OUT/snv_mutated.fa"

echo "== 2/3 converting RNA-editing table to VCF =="
# An A-to-I edit is a transcript-level A>G, so VEP handles it once it is
# given a VCF. Positions come straight from the ANNOVAR table.
python3 - "$RES" > "$OUT/res_sites.vcf" <<'PY'
import csv, sys
print("##fileformat=VCFv4.2")
print("##reference=GRCh38")
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    c = (r.get("Chr") or "").strip()
    if not c or c.upper() == "CHROM":
        continue
    print(f"{c}\t{r['Start']}\t.\t{r['Ref']}\t{r['Alt']}\t.\tPASS\tSOURCE=RES")
PY
if command -v bcftools >/dev/null; then
  bcftools sort "$OUT/res_sites.vcf" -o "$OUT/res_sites.sorted.vcf"
else
  cp "$OUT/res_sites.vcf" "$OUT/res_sites.sorted.vcf"
fi

echo "== 3/3 VEP on RNA-editing sites =="
vep --input_file "$OUT/res_sites.sorted.vcf" --output_file "$OUT/res.txt" \
    --force_overwrite --species homo_sapiens --assembly GRCh38 \
    --cache --offline --dir_cache "$CACHE" --fasta "$GENOME" \
    --symbol --protein --hgvs --pick \
    --plugin ProteinSeqs,"$OUT/res_reference.fa","$OUT/res_mutated.fa"

cat "$OUT"/snv_mutated.fa "$OUT"/res_mutated.fa > "$OUT/vep_all_mutated.fa"

cat <<TXT

VEP output ready. Compare with:

  python3 scripts/05_compare_tools.py \\
    --fasta ours=results/fasta/HCC1395_variant_proteins.uniprot.representative.fasta \\
    --fasta vep=$OUT/vep_all_mutated.fa \\
    --k 9 --outdir results

Expect high 9-mer Jaccard. VEP does not handle fusions or splicing
events, so those genes will show as ours-only; that is a scope
difference, not a disagreement.
TXT

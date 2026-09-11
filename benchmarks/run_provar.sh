#!/usr/bin/env bash
# Run ProVar (ProGenNo/ProHap) on the HCC1395 somatic VCF.
#
# ProVar is the fair comparator for v2p: same input shape (one sample's
# VCF), same scope, same authors as ProHap. ProHap itself needs phased
# population panels and about a terabyte, and comparing tools built for
# different jobs invites the comparison to be dismissed.
#
# This drives src/provar.py per chromosome directly rather than through
# Snakemake. The pipeline's own Snakefile does the same thing; going
# direct avoids installing Snakemake and a conda solve for what is, in
# the end, one Python script per chromosome. The environment it needs is
# small - python, numpy, biopython, gffutils, pandas, pyvcf3, pyarrow -
# and micromamba installs it without root.
#
# Reference: Ensembl 110, the release GENCODE v44 is built from, so both
# tools translate against the same annotation. Anything else would
# confound the tool with its reference.
#
# Usage:  bash benchmarks/run_provar.sh <work_dir> <hcc1395_vcf>
set -euo pipefail

WORK="${1:-/mnt/d/conversion/provar_bench}"
VCF_IN="${2:-/mnt/d/conversion/mnt/data/sSNV_sIndel.vcf.gz}"
PROHAP="$WORK/ProHap"
REF="$WORK/ref"
OUT="$WORK/out"
MM="$WORK/bin/micromamba"
export MAMBA_ROOT_PREFIX="$WORK/mamba"
run() { "$MM" run -n provar "$@"; }

# Spelled out rather than built with seq: seq separates with newlines,
# and this list is interpolated into a nested `bash -c` below, where a
# newline ends the statement and the loop dies with a syntax error.
CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

mkdir -p "$OUT" "$REF/gtf" "$REF/fasta" "$OUT/vcf" "$OUT/tmp" "$OUT/log"

# ---- reference ----------------------------------------------------
if [ ! -s "$REF/gtf/ens110.gtf" ]; then
  echo "[ref] decompressing GTF"
  gunzip -c "$REF/Homo_sapiens.GRCh38.110.chr.gtf.gz" > "$REF/gtf/ens110.gtf"
fi

if [ ! -s "$REF/fasta/total_cdnas_110.fa" ]; then
  echo "[ref] merging cDNA + ncRNA"
  gunzip -c "$REF/Homo_sapiens.GRCh38.ncrna.fa.gz" \
           "$REF/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
    > "$REF/fasta/total_cdnas_110.fa"
fi

# Every protein-coding transcript, so ProVar is not handicapped by a
# narrower transcript set than v2p's --transcript-mode all.
if [ ! -s "$REF/transcripts_110.csv" ]; then
  echo "[ref] building transcript list"
  { echo "chromosome,transcriptID";
    awk -F'\t' '$3=="transcript" && /transcript_biotype "protein_coding"/ {
        match($9, /transcript_id "[^"]+"/);
        id = substr($9, RSTART+15, RLENGTH-16);
        print $1 "," id
      }' "$REF/gtf/ens110.gtf"; } > "$REF/transcripts_110.csv"
fi

# ---- per-chromosome annotation databases ---------------------------
for c in $CHROMS; do
  if [ ! -s "$REF/gtf/ens110_chr${c}.db" ]; then
    echo "[gtf] chromosome $c"
    { grep "^#" "$REF/gtf/ens110.gtf" || true; \
      grep -P "^${c}\t" "$REF/gtf/ens110.gtf" || true; } \
      > "$REF/gtf/ens110_chr${c}.gtf"
    run python3 "$PROHAP/src/parse_gtf.py" \
        -i "$REF/gtf/ens110_chr${c}.gtf" -o "$REF/gtf/ens110_chr${c}.db"
    rm -f "$REF/gtf/ens110_chr${c}.gtf"
  fi
done

# ---- input VCF -----------------------------------------------------
# The truth VCF names contigs chr1..chrX; Ensembl names them 1..X. That
# is a naming convention, not a difference in content, so renaming here
# is not a thumb on the scale.
if [ ! -s "$OUT/vcf/hcc1395.noprefix.vcf" ]; then
  echo "[vcf] stripping chr prefix"
  gunzip -c "$VCF_IN" \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {sub(/^chr/,"",$1); print}' \
    > "$OUT/vcf/hcc1395.noprefix.vcf"
fi

if [ ! -s "$OUT/vcf/ready" ]; then
  echo "[vcf] splitting per chromosome"
  run python3 "$PROHAP/src/fragment_variant_vcf.py" \
      -i "$OUT/vcf/hcc1395.noprefix.vcf" -o "$OUT/vcf/variants"
  touch "$OUT/vcf/ready"
fi

# ---- ProVar --------------------------------------------------------
echo "[provar] translating"
/usr/bin/time -v -o "$OUT/provar_time.txt" bash -c '
  set -e
  for c in '"$CHROMS"'; do
    v="'"$OUT"'/vcf/variants_chr${c}.vcf.gz"
    [ -s "$v" ] || continue
    echo "  chr${c}"
    "'"$MM"'" run -n provar python3 "'"$PROHAP"'/src/provar.py" \
      -i "$v" \
      -db "'"$REF"'/gtf/ens110_chr${c}.db" \
      -transcripts "'"$REF"'/transcripts_110.csv" \
      -cdna "'"$REF"'/fasta/total_cdnas_110.fa" \
      -chr "$c" -acc_prefix "hcc1395_" -af 0.0 -require_start 1 \
      -log "'"$OUT"'/log/provar_chr${c}.log" \
      -tmp_dir "'"$OUT"'/tmp" \
      -output_csv "'"$OUT"'/provar_chr${c}.tsv.gz" \
      -output_fasta "'"$OUT"'/provar_chr${c}.fa"
  done
'

echo "[provar] merging"
cat "$OUT"/provar_chr*.fa > "$OUT/provar_variants.fa"
run python3 - "$OUT" <<'PY'
import glob, gzip, os, sys
out = sys.argv[1]
rows, header = [], None
for f in sorted(glob.glob(os.path.join(out, "provar_chr*.tsv.gz"))):
    # A chromosome with no variants gets ProVar's empty_output(), which
    # writes plain text under the .gz name, so sniff rather than trust
    # the extension.
    with open(f, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    with (gzip.open if gzipped else open)(f, "rt") as fh:
        h = fh.readline()
        header = header or h
        rows.extend(fh.readlines())
with open(os.path.join(out, "provar_variants.tsv"), "w") as fh:
    fh.write(header or "")
    fh.writelines(rows)
print(f"merged {len(rows)} rows")
PY

grep -c '^>' "$OUT/provar_variants.fa" | xargs echo "[provar] sequences:"
echo "[provar] done -> $OUT/provar_variants.{fa,tsv}"

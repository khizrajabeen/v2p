#!/usr/bin/env bash
# Stage 0b - download, index and verify every reference this pipeline needs.
#
#   bash scripts/00_fetch_references.sh ref/
#
# Downloads ~4.5 GB, uses ~18 GB on disk after decompression and indexing.
# Safe to re-run: every step is skipped if its output already exists.
#
# Build choice is not a preference. The RES filename says hg38 and the SEQC2
# truth set is GRCh38, so a GRCh37 reference produces silently wrong codons
# rather than an error. This script asserts the build before you spend an
# hour translating.

set -euo pipefail

REFDIR="${1:-ref}"
GENCODE_RELEASE="${GENCODE_RELEASE:-44}"
mkdir -p "$REFDIR"
cd "$REFDIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "missing: $1"; exit 1; }; }
need curl
need python3

AVAIL_KB=$(df -Pk . | awk 'NR==2{print $4}')
if [ "$AVAIL_KB" -lt 20000000 ]; then
  echo "WARNING: only $((AVAIL_KB/1024/1024)) GB free here; ~18 GB needed." >&2
  echo "Press Ctrl-C to abort, or Enter to continue anyway." >&2
  read -r _
fi

BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
GENOME_GZ="GRCh38.primary_assembly.genome.fa.gz"
GENOME="GRCh38.primary_assembly.genome.fa"
GTF="gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
PROT="gencode.v${GENCODE_RELEASE}.pc_translations.fa.gz"

echo "== 1/5 downloading GENCODE release ${GENCODE_RELEASE} =="
# -C - resumes a partial download; a truncated FASTA is the classic silent failure
[ -f "$GENOME_GZ" ] || curl -fL -C - -O "${BASE}/${GENOME_GZ}"
[ -f "$GTF" ]       || curl -fL -C - -O "${BASE}/${GTF}"
# GENCODE's own translations: a free third opinion for the validation gate
[ -f "$PROT" ]      || curl -fL -C - -O "${BASE}/${PROT}" || \
  echo "  (pc_translations optional; continuing without it)"

echo "== 2/5 decompressing genome =="
# pyfaidx needs plain or bgzip; GENCODE ships plain gzip, which is not
# randomly accessible. Decompress rather than fight it.
[ -f "$GENOME" ] || gunzip -k "$GENOME_GZ"

echo "== 3/5 verifying integrity =="
python3 - "$GTF" <<'PY'
import gzip, sys
n = 0
with gzip.open(sys.argv[1], "rt") as fh:
    for line in fh:
        n += 1
print(f"  GTF readable to the end: {n:,} lines")
PY

echo "== 4/5 indexing genome =="
python3 - "$GENOME" <<'PY'
import sys
try:
    from pyfaidx import Fasta
except ImportError:
    sys.exit("  pyfaidx not installed: pip install pyfaidx")
fa = Fasta(sys.argv[1], sequence_always_upper=True)
names = list(fa.keys())
print(f"  contigs: {len(names)}   first: {names[:3]}")
print(f"  naming: {'UCSC (chr-prefixed)' if names[0].startswith('chr') else 'Ensembl'}")
PY

echo "== 5/5 build assertion =="
python3 - "$GENOME" <<'PY'
import sys
from pyfaidx import Fasta
fa = Fasta(sys.argv[1], sequence_always_upper=True)
EXPECT = {"chr1": 248956422, "chr2": 242193529, "chr17": 83257441}
GRCH37 = {"chr1": 249250621}
bad = []
for c, want in EXPECT.items():
    key = c if c in fa else c.replace("chr", "")
    if key not in fa:
        bad.append(f"{c}: absent")
        continue
    got = len(fa[key])
    if got != want:
        hint = " (this is GRCh37/hg19)" if got == GRCH37.get(c) else ""
        bad.append(f"{c}: {got:,} != {want:,}{hint}")
    else:
        print(f"  {c} = {got:,} bp  OK")
if bad:
    sys.exit("  BUILD MISMATCH:\n    " + "\n    ".join(bad))
print("  GRCh38/hg38 confirmed")
PY

if [ ! -f uniprot_human_SP.fasta ]; then
  cat <<'EOF'

NOTE: uniprot_human_SP.fasta is not here.
Copy the file Zhiyin supplied into this directory:

    cp Human_Homo_sapiens_uniprot_SP_UP000005640_N20432_20260304.fasta \
       ref/uniprot_human_SP.fasta

It is the format template AND the validation ground truth - the pipeline
uses it for both.
EOF
fi

cat <<EOF

References ready in $PWD

  --genome  $PWD/${GENOME}
  --gtf     $PWD/${GTF}
  --uniprot $PWD/uniprot_human_SP.fasta

Sanity-check them together with:
  python3 ../scripts/00_inspect_inputs.py ../data $PWD
EOF

#!/usr/bin/env bash
# Stage 0b - download, index and verify a reference for this pipeline.
#
#   bash scripts/00_fetch_references.sh ref/                     # GRCh38
#   bash scripts/00_fetch_references.sh ref37/  --assembly GRCh37
#   bash scripts/00_fetch_references.sh reft2t/ --assembly T2T
#
# Downloads ~4.5 GB and uses ~18 GB on disk after decompression and
# indexing. Safe to re-run: every step is skipped if its output exists.
#
# The assembly is asserted by contig length before you spend an hour
# translating. Using GRCh37 against GRCh38 calls produces silently wrong
# codons rather than an error, which is why that check exists and runs
# before anything else is trusted.
#
# Assemblies
#   GRCh38  GENCODE primary assembly and annotation. The default.
#   GRCh37  GENCODE's own GRCh37 mapping of the same release, so the
#           annotation release matches across builds and only the
#           coordinates differ.
#   T2T     CHM13v2.0 from UCSC (hs1), chr-named so it lines up with
#           chr-prefixed calls, annotated with NCBI RefSeq lifted to T2T.
#
# Every URL below was checked to resolve before being written down.

set -euo pipefail

REFDIR=""
ASSEMBLY="GRCh38"
while [ $# -gt 0 ]; do
  case "$1" in
    --assembly)   ASSEMBLY="${2:-}"; shift 2 ;;
    --assembly=*) ASSEMBLY="${1#*=}"; shift ;;
    -h|--help)    sed -n '2,25p' "$0"; exit 0 ;;
    *)            REFDIR="$1"; shift ;;
  esac
done
REFDIR="${REFDIR:-ref}"
GENCODE_RELEASE="${GENCODE_RELEASE:-44}"

case "$ASSEMBLY" in
  GRCh38|grch38|hg38)      ASSEMBLY=GRCh38 ;;
  GRCh37|grch37|hg19)      ASSEMBLY=GRCh37 ;;
  T2T|t2t|CHM13|chm13|hs1) ASSEMBLY=T2T ;;
  *) echo "unknown assembly: $ASSEMBLY (choose GRCh38, GRCh37 or T2T)" >&2
     exit 2 ;;
esac

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

GENCODE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips"
PROT_URL=""

case "$ASSEMBLY" in
  GRCh38)
    GENOME_URL="${GENCODE}/GRCh38.primary_assembly.genome.fa.gz"
    GTF_URL="${GENCODE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
    PROT_URL="${GENCODE}/gencode.v${GENCODE_RELEASE}.pc_translations.fa.gz"
    GENOME="GRCh38.primary_assembly.genome.fa"
    ;;
  GRCh37)
    M="${GENCODE}/GRCh37_mapping"
    GENOME_URL="${M}/GRCh37.primary_assembly.genome.fa.gz"
    GTF_URL="${M}/gencode.v${GENCODE_RELEASE}lift37.annotation.gtf.gz"
    PROT_URL="${M}/gencode.v${GENCODE_RELEASE}lift37.pc_translations.fa.gz"
    GENOME="GRCh37.primary_assembly.genome.fa"
    ;;
  T2T)
    # UCSC's hs1 rather than NCBI's release: NCBI names contigs by RefSeq
    # accession (NC_060925.1), which will not match a chr-prefixed VCF,
    # and a contig-naming mismatch fails the whole run.
    GENOME_URL="${UCSC}/hs1.fa.gz"
    GTF_URL="${UCSC}/genes/hs1.ncbiRefSeq.gtf.gz"
    GENOME="hs1.fa"
    ;;
esac

GENOME_GZ="$(basename "$GENOME_URL")"
GTF="$(basename "$GTF_URL")"
PROT="${PROT_URL:+$(basename "$PROT_URL")}"

echo "== 1/5 downloading ${ASSEMBLY} (GENCODE release ${GENCODE_RELEASE}) =="
# -C - resumes a partial download; a truncated FASTA is the classic silent
# failure this script exists to prevent.
[ -f "$GENOME" ] || [ -f "$GENOME_GZ" ] || curl -fL -C - -O "$GENOME_URL"
[ -f "$GTF" ] || curl -fL -C - -O "$GTF_URL"
if [ -n "$PROT_URL" ]; then
  [ -f "$PROT" ] || curl -fL -C - -O "$PROT_URL" || \
    echo "  (translations optional; continuing without them)"
fi

echo "== 2/5 decompressing genome =="
# pyfaidx needs plain or bgzip; these ship plain gzip, which is not
# randomly accessible. Decompress rather than fight it.
[ -f "$GENOME" ] || gunzip -k "$GENOME_GZ"

echo "== 3/5 verifying integrity =="
python3 - "$GTF" <<'PY'
import gzip, sys
n = 0
with gzip.open(sys.argv[1], "rt") as fh:
    for _ in fh:
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
python3 - "$GENOME" "$ASSEMBLY" <<'PY'
import sys
from pyfaidx import Fasta

# Contig lengths taken from each assembly's own chrom.sizes rather than
# from memory. They are the cheapest unambiguous discriminator there is.
LENGTHS = {
    "GRCh38": {"chr1": 248956422, "chr2": 242193529, "chr17": 83257441},
    "GRCh37": {"chr1": 249250621, "chr2": 243199373},
    "T2T":    {"chr1": 248387328, "chr2": 242696752, "chr3": 201105948},
}
OTHERS = {a: v["chr1"] for a, v in LENGTHS.items()}

fa_path, assembly = sys.argv[1], sys.argv[2]
fa = Fasta(fa_path, sequence_always_upper=True)
expect = LENGTHS[assembly]

bad = []
for c, want in expect.items():
    key = c if c in fa else c.replace("chr", "")
    if key not in fa:
        bad.append(f"{c}: absent")
        continue
    got = len(fa[key])
    if got == want:
        print(f"  {c} = {got:,} bp  OK")
        continue
    # Naming the assembly it actually is turns "wrong length" into an
    # instruction; the usual cause is fetching one build over another.
    hint = ""
    if c == "chr1":
        match = [a for a, ln in OTHERS.items() if ln == got]
        if match:
            hint = f" (this is {match[0]})"
    bad.append(f"{c}: {got:,} != {want:,}{hint}")

if bad:
    sys.exit(f"  BUILD MISMATCH - expected {assembly}:\n    "
             + "\n    ".join(bad))
print(f"  {assembly} confirmed")
PY

if ! ls ./*uniprot*.fasta ./*sprot*.fasta >/dev/null 2>&1; then
  cat <<'EOF'

NOTE: no reference proteome is in this directory.
Download the SwissProt set for your organism from UniProt and put it here,
for example:

    curl -L -o uniprot_human_SP.fasta \
      'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'

It serves two purposes: the format template the output matches, and the
proteome the validation stage compares against.
EOF
fi

if [ -z "$PROT_URL" ]; then
  cat <<'EOF'

NOTE: this assembly ships no "annotation source's own translations" file.
That file drives the strongest validation gate - our translation of a
transcript against the annotation's own translation of the same
transcript. Without it, validation falls back to comparing against
UniProt canonical, which disagrees on start codons often enough to fail
the default threshold for reasons that are not errors. Either lower the
gate with --min-agreement, or supply a translations FASTA with
--translations.
EOF
fi

cat <<EOF

Reference ready in $PWD  (${ASSEMBLY})

  --genome      $PWD/${GENOME}
  --annotation  $PWD/${GTF}
${PROT:+  --translations $PWD/${PROT}}
  --proteome    <your SwissProt FASTA>

Or pass --ref $PWD and let detection work it out.
EOF

#!/usr/bin/env bash
# Build the HCC1395 BENCHMARK protein database — one command, no arguments.
#
#   bash benchmark.sh
#
# Run from the project root, with references in ref/ and inputs in data/.
#
# Benchmark means maximum coverage, not a minimal search database:
#   every variant class, synonymous and UTR variants kept, NMD-flagged
#   entries kept, every coding transcript, and the complete SwissProt
#   proteome present verbatim under its own accessions.

set -euo pipefail
cd "$(dirname "$0")"

GENOME=ref/GRCh38.primary_assembly.genome.fa
GTF=ref/gencode.v44.annotation.gtf.gz
UNIPROT=ref/uniprot_human_SP.fasta
GENCODE_PROT=ref/gencode.v44.pc_translations.fa.gz
OUT=benchmark
NAME=HCC1395_benchmark_proteome

for f in "$GENOME" "$GTF" "$UNIPROT"; do
  [ -f "$f" ] || { echo "missing reference: $f"; exit 1; }
done

echo "=== 1/6  parse inputs ==="
python3 src/v2p/stages/01_parse_inputs.py \
    --vcf    data/sSNV_sIndel.vcf.gz \
    --res    data/HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt \
    --fusion data/HCC1395_high_confidence_Fusion_genes_all.csv \
    --as-lr  data/HCC1395_high_confidence_AS-LR_v1.csv \
    --outdir results >/dev/null 2>&1

echo "=== 2/6  input QC ==="
python3 src/v2p/stages/03_qc_report.py \
    --manifest results/tables/unified_variant_manifest.tsv \
    --outdir results >/dev/null 2>&1

echo "=== 3/6  translate (all transcripts, keeping everything) ==="
python3 src/v2p/stages/02_build_protein_fasta.py \
    --manifest results/tables/unified_variant_manifest.tsv \
    --genome "$GENOME" --gtf "$GTF" --uniprot "$UNIPROT" \
    --header-style uniprot \
    --transcript-mode all \
    --keep-synonymous \
    --include-reference \
    --emit-disposition results/tables/disposition.benchmark.tsv \
    --emit-isoform-gtf results/tables/as_isoforms.benchmark.gtf \
    --outdir results 2>&1 | grep -E "wrote .* sequences" || true

BUILT=results/fasta/HCC1395_variant_proteins.uniprot.all.fasta

echo "=== 4/6  validate ==="
VALARGS=(--uniprot "$UNIPROT"
         --recoding results/tables/res_recoding_sites.tsv
         --protein-fasta "$BUILT" --outdir results)
[ -f "$GENCODE_PROT" ] && VALARGS+=(--gencode-translations "$GENCODE_PROT")
python3 src/v2p/stages/04_validate_uniprot.py "${VALARGS[@]}" >/dev/null 2>&1 || {
  echo "VALIDATION FAILED — see results/qc/uniprot_validation.md"; exit 1; }
grep -E "Exact-agreement|Same-protein" results/qc/uniprot_validation.md || true

python3 src/v2p/stages/08_recovery_report.py \
    --manifest results/tables/unified_variant_manifest.tsv \
    --disposition results/tables/disposition.benchmark.tsv \
    --outdir results >/dev/null 2>&1

echo "=== 5/6  package ==="
rm -rf "$OUT"
python3 src/v2p/stages/06_package_release.py \
    --fasta "$BUILT" --uniprot "$UNIPROT" \
    --name "$NAME" --outdir "$OUT" \
    --decoy pseudo_reverse --split-by-class --append-reference \
    --extra-file results/qc >/dev/null 2>&1

echo "=== 6/6  complete the reference proteome and rebuild decoys ==="
# The packager skips a UniProt entry whose sequence already appears among
# the variant records. That is right for a minimal search database and
# wrong for a benchmark: those proteins would then exist only under a
# vr| header, so a search result could not be mapped back to a standard
# accession. Append them, then regenerate the decoys over the full set.
python3 - "$OUT" "$NAME" "$UNIPROT" <<'PY'
import re, sys
from pathlib import Path

OUT, NAME, UNIPROT = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
tgt = OUT / f"{NAME}.target.fasta"

def read(p):
    h, buf = None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h: yield h, "".join(buf)
            h, buf = line, []
        elif h: buf.append(line.strip())
    if h: yield h, "".join(buf)

cur  = list(read(tgt))
have = {h.split()[0] for h, _ in cur}
add  = [(h, s) for h, s in read(UNIPROT) if h.split()[0] not in have]
allr = cur + add
print(f"  appended {len(add)} SwissProt entries the packager had skipped")

with open(tgt, "w") as fh:
    for h, s in allr: fh.write(h + "\n" + s + "\n")

with open(OUT / "by_class" / f"{NAME}.REFERENCE_PROTEOME.fasta", "w") as fh:
    for h, s in allr:
        if h.startswith(">sp|"): fh.write(h + "\n" + s + "\n")

TRYP = re.compile(r"(?<=[KR])(?!P)")
def pseudo_reverse(seq):
    out = []
    for f in TRYP.split(seq):
        if not f: continue
        out.append(f[:-1][::-1] + f[-1] if f[-1] in "KR" else f[::-1])
    return "".join(out)

with open(OUT / f"{NAME}.target_decoy.fasta", "w") as fh:
    for h, s in allr: fh.write(h + "\n" + s + "\n")
    for h, s in allr:
        first, *rest = h[1:].split(" ", 1)
        fh.write(f">DECOY_{first}" + (" " + rest[0] if rest else "")
                 + " DECOY=pseudo_reverse\n")
        fh.write(pseudo_reverse(s) + "\n")
PY

mkdir -p "$OUT/tables"
cp results/tables/disposition.benchmark.tsv "$OUT/tables/" 2>/dev/null || true
cp results/tables/recovery_by_class.tsv     "$OUT/tables/" 2>/dev/null || true

cat > "$OUT/README.md" <<'EOF'
# HCC1395 benchmark protein database

Complete representation of every variant that yields a protein. Nothing is
filtered for novelty, decay or redundancy.

Start with `HCC1395_benchmark_proteome.target.fasta` — all variant types in
one file. Every entry carries `VT=<type>` in its header.

- `*.target_decoy.fasta` — the same plus matched pseudo-reverse decoys
- `*.entries.tsv` — one row per sequence, including the sequence itself
- `*.summary_by_class.tsv` — composition per variant type
- `by_class/` — the same entries split by type; `INDEX.md` explains each
- `tables/disposition.benchmark.tsv` — one row per input variant, so the
  full 42,369-variant-to-sequence mapping is recoverable
- `METHODS.md` — construction, validation and caveats
- `qc/recovery_report.md` — the fate of every input variant
- `qc/uniprot_validation.md` — agreement with GENCODE and UniProt
- `MANIFEST.txt` — SHA-256 of every file

Verify integrity:
  awk 'NF>=3{print $1"  "$2}' MANIFEST.txt | sha256sum -c -

Included on purpose: synonymous and UTR variants (protein identical to the
reference), NMD-flagged truncations, every coding transcript rather than one
per gene, and the complete SwissProt proteome under its own accessions.

Where several variants give the same amino-acid sequence they appear as one
entry, with the other classes listed in `ALSO=`. A FASTA cannot hold the same
sequence twice in a way a search engine can distinguish; the per-variant
detail is in `tables/disposition.benchmark.tsv`.
EOF

( cd "$OUT" && { echo "# MANIFEST — $NAME"; echo
    find . -type f ! -name MANIFEST.txt -printf '%P\n' | sort |
    while read -r f; do
      echo "$(sha256sum "$f" | cut -d' ' -f1)  $f  ($(stat -c%s "$f") bytes)"
    done; } > MANIFEST.txt )

echo
echo "done -> $OUT/"
printf "  target sequences   %s\n" "$(grep -c '^>' "$OUT/$NAME.target.fasta")"
printf "  sp| accessions     %s\n" "$(grep -c '^>sp|' "$OUT/$NAME.target.fasta")"
printf "  decoys             %s\n" "$(grep -c '^>DECOY_' "$OUT/$NAME.target_decoy.fasta")"
printf "  by_class files     %s\n" "$(ls "$OUT/by_class" | wc -l)"

#!/usr/bin/env bash
# Rename the built output into descriptive, variant-type-based filenames
# and convert every .md to plain .txt.
#
#   bash finalize.sh
#
# Run from the project root after benchmark.sh. Reads benchmark/ and
# writes HCC1395_variant_proteins/ — the folder to send.

set -euo pipefail
cd "$(dirname "$0")"

SRC=benchmark
OLD=HCC1395_benchmark_proteome
DST=HCC1395_variant_proteins
PFX=HCC1395

[ -d "$SRC" ] || { echo "no $SRC/ — run benchmark.sh first"; exit 1; }
rm -rf "$DST"; mkdir -p "$DST/by_variant_type" "$DST/tables" "$DST/logs"

# variant class -> readable file name
declare -A LABEL=(
  [SNV]=SNV
  [INDEL]=InDel
  [RNA_EDITING]=RNA_editing
  [FUSION]=gene_fusion
  [AS_SE]=splicing_skipped_exon
  [AS_RI]=splicing_retained_intron
  [AS_MX]=splicing_mutually_exclusive_exons
  [AS_A3]=splicing_alt_3prime_splice_site
  [AS_A5]=splicing_alt_5prime_splice_site
  [AS_AF]=splicing_alt_first_exon
  [AS_AL]=splicing_alt_last_exon
  [REFERENCE]=wildtype_counterparts
  [REFERENCE_PROTEOME]=reference_proteome_SwissProt
)

echo "== combined =="
cp "$SRC/$OLD.target.fasta"       "$DST/${PFX}_all_variant_types.fasta"
cp "$SRC/$OLD.target_decoy.fasta" "$DST/${PFX}_all_variant_types.with_decoys.fasta"
cp "$SRC/$OLD.entries.tsv"        "$DST/${PFX}_entry_details.tsv"
cp "$SRC/$OLD.summary_by_class.tsv" "$DST/${PFX}_summary_by_variant_type.tsv"

echo "== per variant type =="
for f in "$SRC"/by_variant_type/*.fasta "$SRC"/by_class/*.fasta; do
  [ -e "$f" ] || continue
  cls=$(basename "$f" .fasta); cls=${cls#"$OLD."}
  out=${LABEL[$cls]:-$cls}
  cp "$f" "$DST/by_variant_type/${PFX}_${out}.fasta"
  printf "  %-34s %s\n" "${PFX}_${out}.fasta" "$(grep -c '^>' "$f")"
done

echo "== tables and logs =="
cp "$SRC"/tables/*.tsv "$DST/tables/" 2>/dev/null || true
cp "$SRC"/logs/*        "$DST/logs/"   2>/dev/null || true

# markdown -> plain text: strip heading hashes, bold markers, backticks,
# table pipes and bullet dashes, so the files read cleanly in Notepad.
md2txt () {
  sed -e 's/^#\{1,6\} *//' \
      -e 's/\*\*//g' -e 's/`//g' \
      -e 's/^| *//' -e 's/ *|$//' -e 's/ *| */  /g' \
      -e '/^[- :]\{4,\}$/d' \
      -e 's/^- /  * /' "$1" > "$2"
}

echo "== documentation as .txt =="
[ -f "$SRC/METHODS.md" ]                 && md2txt "$SRC/METHODS.md"                 "$DST/METHODS.txt"
[ -f "$SRC/qc/recovery_report.md" ]      && md2txt "$SRC/qc/recovery_report.md"      "$DST/RECOVERY_REPORT.txt"
[ -f "$SRC/qc/uniprot_validation.md" ]   && md2txt "$SRC/qc/uniprot_validation.md"   "$DST/VALIDATION_REPORT.txt"
[ -f "$SRC/qc/summary.md" ]              && md2txt "$SRC/qc/summary.md"              "$DST/INPUT_SUMMARY.txt"

N_ALL=$(grep -c '^>' "$DST/${PFX}_all_variant_types.fasta")
N_SP=$(grep -c '^>sp|' "$DST/${PFX}_all_variant_types.fasta")
N_DEC=$(grep -c '^>DECOY_' "$DST/${PFX}_all_variant_types.with_decoys.fasta")

cat > "$DST/README.txt" <<EOF
HCC1395 VARIANT PROTEIN SEQUENCES
=================================

Amino-acid sequences for every variant in the HCC1395 high-confidence
package: somatic SNVs and InDels, A-to-I RNA editing sites, gene fusions,
and long-read alternative splicing events.

START HERE
----------
${PFX}_all_variant_types.fasta
    $N_ALL sequences, all variant types in one file. Every entry carries
    VT=<variant type> in its header, alongside the consequence, transcript,
    protein change, genomic locus and NMD prediction.

    Header format follows the SwissProt template, sequences unwrapped:
      >vr|<id>|<GENE>_HUMAN_<TYPE> <description> OS=Homo sapiens OX=9606
       GN=<gene> VT=<type> CSQ=<consequence> TX=<transcript> PC=<change> ...

    Reference entries keep the standard sp| prefix; variant entries use
    vr|, so the two are separable on the database field alone.

OTHER FILES
-----------
${PFX}_all_variant_types.with_decoys.fasta
    The same plus $N_DEC pseudo-reverse decoys, for false-discovery-rate
    control in a database search.

${PFX}_entry_details.tsv
    One row per sequence: variant type, gene, transcript, protein change,
    locus, NMD flag, tryptic peptide counts, and the sequence itself.

${PFX}_summary_by_variant_type.tsv
    Counts per variant type.

by_variant_type/
    The same sequences split into one file per variant type. See
    FILE_INDEX.txt for what each contains.

tables/
    disposition — one row per input variant, recording whether it produced
    a protein and if not, why. The full 42,369-variant mapping.
    recovery_by_class — recovery counts per variant type.

METHODS.txt            how the sequences were built and validated
RECOVERY_REPORT.txt    the fate of every input variant
VALIDATION_REPORT.txt  agreement with GENCODE and UniProt
INPUT_SUMMARY.txt      composition of the input variant set
logs/                  per-stage provenance, including input checksums
MANIFEST.txt           SHA-256 of every file here

SCOPE
-----
Nothing is filtered for novelty or redundancy. Included on purpose:
synonymous and UTR variants whose protein equals the reference, entries
flagged as likely nonsense-mediated-decay targets, every coding transcript
rather than one per gene, and the complete SwissProt reference proteome
($N_SP sp| entries, comprising the 20,432 UniProt proteins plus the
wild-type counterparts of variant-carrying transcripts).

Where several variants give the same amino-acid sequence they appear as a
single entry, with the other types listed in ALSO=. A FASTA cannot hold the
same sequence twice in a way a search engine can distinguish; the
per-variant detail is kept in tables/disposition.

VERIFY
------
  awk 'NF>=3{print \$1"  "\$2}' MANIFEST.txt | sha256sum -c -
EOF

{
  echo "VARIANT TYPE FILES"
  echo "=================="
  echo
  echo "Each file is a subset of ${PFX}_all_variant_types.fasta with the same"
  echo "headers. Together they contain every sequence, with no overlap."
  echo
  for cls in SNV INDEL RNA_EDITING FUSION AS_SE AS_RI AS_MX AS_A3 AS_A5 \
             AS_AF AS_AL REFERENCE REFERENCE_PROTEOME; do
    out=${LABEL[$cls]}
    fp="$DST/by_variant_type/${PFX}_${out}.fasta"
    [ -f "$fp" ] || continue
    printf "%s_%s.fasta\n" "$PFX" "$out"
    printf "    %s sequences.  Header field: VT=%s\n" "$(grep -c '^>' "$fp")" "$cls"
    case $cls in
      SNV) echo "    Somatic single-nucleotide variants: one amino-acid substitution,";
           echo "    a premature stop, or a lost stop.";;
      INDEL) echo "    Somatic insertions and deletions. In-frame changes alter residues";
             echo "    locally; frameshifts replace the entire C-terminus.";;
      RNA_EDITING) echo "    A-to-I editing. ADAR converts adenosine to inosine, read as";
                   echo "    guanosine. The DNA is unchanged, so these proteoforms exist";
                   echo "    only at the RNA and protein level.";;
      FUSION) echo "    Chimeric proteins: the 5' partner up to the breakpoint joined to";
              echo "    the 3' partner after it. The header records whether the junction";
              echo "    keeps the codon boundary and whether the 3' partner stays in its";
              echo "    own reading frame.";;
      AS_SE) echo "    Skipped exon: a cassette exon included or excluded.";;
      AS_RI) echo "    Retained intron: an intron kept in the mature transcript, usually";
             echo "    introducing a premature stop. Check the NMD flag.";;
      AS_MX) echo "    Mutually exclusive exons: exactly one of two alternatives used.";;
      AS_A3) echo "    Alternative 3' splice site: a different acceptor shifts the exon start.";;
      AS_A5) echo "    Alternative 5' splice site: a different donor shifts the exon end.";;
      AS_AF) echo "    Alternative first exon: a different start codon and N-terminus.";
             echo "    The largest class in this dataset.";;
      AS_AL) echo "    Alternative last exon: a different C-terminus.";;
      REFERENCE) echo "    Wild-type protein of each transcript carrying a variant. A variant";
                 echo "    peptide is only interpretable next to its wild-type counterpart.";;
      REFERENCE_PROTEOME) echo "    The supplied UniProt SwissProt entries, verbatim and unmodified.";;
    esac
    echo
  done
} > "$DST/FILE_INDEX.txt"

( cd "$DST" && { echo "MANIFEST - HCC1395 variant protein sequences"; echo
    find . -type f ! -name MANIFEST.txt -printf '%P\n' | sort |
    while read -r f; do
      echo "$(sha256sum "$f" | cut -d' ' -f1)  $f  ($(stat -c%s "$f") bytes)"
    done; } > MANIFEST.txt )

echo
echo "done -> $DST/"
find "$DST" -maxdepth 1 -type f -printf '  %-46f %s bytes\n' | sort

# HCC1395 variant → amino-acid sequence conversion

Converts the HCC1395 high-confidence variant package into a single protein
FASTA in which every entry carries its variant-type annotation.

Handles four evidence types:

| input file | variant type | records received |
|---|---|---:|
| `high-confidence_sSNV_sIndel_v1...vcf` | somatic SNV / InDel | **not yet received** |
| `HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt` | RNA editing (A-to-I) | 8,093 |
| `HCC1395_high_confidence_Fusion_genes_all.csv` | gene fusion | 28 |
| `HCC1395_high_confidence_AS-LR_v1.csv` | alternative splicing | 1,027 |

---

## Status

| stage | needs reference? | state |
|---|---|---|
| 0 `00_fetch_references.sh` | — | ready to run |
| 1 `01_parse_inputs.py` | no | **already run**, 9,148 records in `results/tables/` |
| 2 `02_build_protein_fasta.py` | **yes** (hg38 + GENCODE) | code complete, 28/28 tests pass, integration-tested end to end |
| 3 `03_qc_report.py` | no | **already run**, `results/qc/summary.md` |
| 4 `04_validate_uniprot.py` | UniProt only | check A **already run** (19/25 sites anchored); check B pending stage 2 |
| 5 `05_compare_tools.py` | tool outputs | ready; see `docs/PLAN.md` |

Stage 2 is the only step still to execute and it needs a ~4 GB reference
download, so it must run on your machine rather than in this session.

---

## Quick start

```bash
pip install -r requirements.txt
python tests/test_pipeline.py          # 28 assertions, ~2 s, no reference needed
bash scripts/00_fetch_references.sh ref/
make all GENOME=ref/GRCh38.primary_assembly.genome.fa \
         GTF=ref/gencode.v44.annotation.gtf.gz \
         VCF=data/high-confidence_sSNV_sIndel_v1.sort.final.combined.sort.vcf
```

Or step by step:

```bash
python scripts/01_parse_inputs.py \
    --vcf    data/high-confidence_sSNV_sIndel_v1.sort.final.combined.sort.vcf \
    --res    data/HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt \
    --fusion data/HCC1395_high_confidence_Fusion_genes_all.csv \
    --as-lr  data/HCC1395_high_confidence_AS-LR_v1.csv \
    --outdir results

python scripts/03_qc_report.py --manifest results/tables/unified_variant_manifest.tsv

python scripts/02_build_protein_fasta.py \
    --manifest results/tables/unified_variant_manifest.tsv \
    --genome   ref/GRCh38.primary_assembly.genome.fa \
    --gtf      ref/gencode.v44.annotation.gtf.gz \
    --header-style peff \
    --include-reference \
    --outdir   results
```

---

## Method

### Small variants (SNV / MNV / InDel) and RNA editing

The variant is applied to the **mature transcript** sequence, not to genomic
sequence. One code path therefore serves both DNA variants and RNA editing:
an A-to-I edit is indistinguishable from a genomic A>G at the transcript
level, and a minus-strand gene's genomic `T>C` becomes a transcript-level
`A>G` automatically once exon blocks are reverse-complemented. (8,093 of
8,093 sites in the RES file are A>G or T>C, exactly as A-to-I editing
predicts — see `results/qc/summary.md`.)

Steps: locate the variant on the transcript → verify the observed reference
base matches the call (mismatches are **flagged, never silently accepted**) →
splice in the alternative allele → translate from the annotated start codon
to the first stop. Frameshifts translate past the annotated stop into the
3′UTR, which is where the novel peptide lives.

Consequences assigned: `synonymous`, `missense`, `stop_gained`, `stop_lost`,
`start_lost`, `inframe_insertion`, `inframe_deletion`, `frameshift`.
Variants that are intronic or that span a splice junction on a given
transcript yield no protein for that transcript.

### Fusions

The table gives genomic breakpoints without partner strand or transcript
context, so the builder picks a representative coding transcript per partner,
takes the 5′ partner's transcript up to the breakpoint and the 3′ partner's
from the breakpoint onward. A breakpoint inside an intron **snaps to the
flanking exon boundary**, which is the junction that would actually be
spliced. Both partner orientations are attempted, because fusion tables do
not consistently order partners 5′→3′; downstream you can keep whichever
orientation produces a viable ORF.

Two frame flags are reported separately, because they are not the same thing:

- `junction_on_codon_boundary` — the 5′ partner's reading frame is intact at
  the junction;
- `3p_native_frame` / `3p_frameshifted` — whether the 3′ partner continues in
  **its own** native frame. This is the flag that matters for a fusion
  neoantigen; only when both hold is the record labelled `in_frame_fusion`.

### Alternative splicing

Two strategies, in order:

1. **Junction matching (preferred).** Each SUPPA2 event id defines the intron
   junctions separating the two alternative forms. Annotated coding
   transcripts of that gene containing the required junctions are selected
   and translated as-is. Nothing is invented — the sequences are real
   GENCODE isoforms, which keeps spurious peptides out of the search space.
2. **De novo construction (fallback, SE/RI/A3/A5 only).** If no annotated
   transcript carries the junction — the interesting, tumour-specific case —
   the representative transcript's exon chain is edited to realise the event
   and translated from the first viable ATG. These records are flagged
   `constructed_isoform` and carry `constructed=1` so you can filter or
   weight them separately.

SUPPA2 coordinate grammar implemented (validated against the real event ids
in your file):

```
SE  <e1>-<s2> : <e2>-<s3>
A5  <e1>-<s3> : <e2>-<s3>
A3  <e1>-<s2> : <e1>-<s3>
RI  <s1> : <e1>-<s2> : <e2>
MX  <e1>-<s2> : <e2>-<s4> : <e1>-<s3> : <e3>-<s4>
AF  <s1> : <e1>-<s3> : <s2> : <e2>-<s3>
AL  <e1>-<s2> : <e2> : <e1>-<s3> : <e3>
```

Group arity is checked against the grammar; a mismatch raises a warning
rather than being misread.

---

## Output format

Four interchangeable header styles (`--header-style`). **`uniprot` is the
default choice for this project**, since it reproduces the grammar and the
single-line layout of the supplied `Human_Homo_sapiens_uniprot_SP_...fasta`
template exactly. All four carry the variant type. Full field reference in [`docs/FORMAT_SPEC.md`](docs/FORMAT_SPEC.md).

**`peff`** — PSI Extended FASTA Format, the ratified HUPO-PSI standard for
encoding sequence variants in FASTA. Supported by Comet, neXtProt, UniProtKB:

```
>HCC1395:GRIA2_NM_000826_chr4_157336723_A_G \DbUniqueId=... \PName=GRIA2 RNA editing site (A-to-I)
 \GName=GRIA2 \TaxName=Homo sapiens \NcbiTaxId=9606 \Length=883 \VariantSimple=(607|R)
 \VariantType=RNA_EDITING \Consequence=missense \TranscriptId=NM_000826
 \ProteinChange=p.Q607R \GenomicLocus=chr4:157336723A>G \NovelSpan=607-607
 \EvidenceSource=RES_ANNOVAR
```

**`descriptive`** — pipe/key-value, safest default with arbitrary search engines:

```
>HCC1395|GRIA2_NM_000826_chr4_157336723_A_G VT=RNA_EDITING CSQ=missense GN=GRIA2 ...
```

**`pvac`** — compact neoantigen-pipeline style: `>MT.GRIA2.NM_000826.RNA_EDITING.missense.Q607R`

A companion `results/tables/protein_records.tsv` gives one row per FASTA
entry for joins and QC.

`--include-reference` additionally emits the unmodified reference protein of
each affected transcript. Recommended for MS database search: without the
wild-type counterpart you cannot tell a genuine variant peptide from a
mis-assigned one.

---

## Reproducibility

Every stage writes two files to `logs/`:

- `<stage>.<run_id>.log` — timestamped operation log
- `<stage>.<run_id>.provenance.json` — command line, working directory,
  Python version, platform, git commit, full `pip freeze`, all parameters,
  and **SHA-256 of every input and output**

Determinism: no randomness anywhere in the pipeline. Transcript selection
ties break on `MANE_Select` → `Ensembl_canonical` → `basic` tag → longest
CDS → longest transcript → lexicographic transcript id, so the same inputs
always give byte-identical output. The genetic code is hard-coded in
`src/v2p/seqops.py` rather than imported, so a library update cannot silently
change a translation.

Input checksums recorded for this run:

```
5b8245c5...  HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt
232b7192...  HCC1395_high_confidence_Fusion_genes_all.csv
064ef665...  HCC1395_high_confidence_AS-LR_v1.csv
```

---

## Testing

`python tests/test_pipeline.py` — 28 assertions, no reference download needed.

`tests/make_fixture.py` builds a miniature genome with three genes (one `+`
strand, one `−` strand, one fusion partner). The transcript is constructed
first and then *placed* into the contigs, so the expected protein for every
test case is known by construction rather than by trusting the code. Covered:
plus- and minus-strand missense, minus-strand A-to-I recoding, synonymous,
stop-gained, 1 bp frameshift, in-frame 3 bp deletion and insertion, intronic
null result, reference-mismatch flagging, fusion junction indexing, intronic
breakpoint snapping, 3′-partner frame detection, SUPPA2 grammar parsing, and
FASTA round-trip in all three header styles.

---

## Limitations — read before using the output

1. **Isoform choice.** Default `--transcript-mode representative` gives one
   protein per gene. `--transcript-mode all` covers every coding transcript;
   it multiplies database size by roughly 4–6× and inflates the FDR of a
   downstream MS search. Choose deliberately.
2. **Alternative splicing is transcript-level, not proteoform-level.** AF
   events (624 of 1,027, the majority) change the N-terminus and therefore
   the start codon; where no annotated transcript matches, the start codon is
   inferred as the first viable ATG and flagged
   `start_codon_inferred_first_ATG`. Treat those with caution.
3. **RNA-editing haplotypes are not phased.** Each site is applied
   independently. Two edits in one codon would need co-occurrence data that
   is not in the input file.
4. **No NMD prediction.** A frameshift or stop-gain >50 nt upstream of the
   last exon-exon junction is likely degraded and may never yield protein.
   The records are emitted anyway and flagged; filter downstream if you want
   NMD-escaping products only.
5. **Selenoproteins truncate at UGA.** Table 1 is applied without
   recoding, so any selenoprotein is cut short at the first in-frame UGA.
6. **Reference build must be GRCh38.** The stage-0 script asserts
   `chr1 == 248,956,422 bp` before you spend an hour translating; a GRCh37
   reference produces silently wrong codons rather than an error.
7. **Fusion transcripts are representative reconstructions.** They are built
   from breakpoints plus annotation, not from assembled long reads. If the
   long-read assemblies behind `Type=LR` calls exist, translating those
   directly is more accurate.

---

## Layout

```
src/v2p/
  provenance.py     run logging, checksums, environment capture
  seqops.py         genetic code, translation, ORF and HGVS utilities
  annotation.py     GTF transcript models, genome access, coordinate mapping
  fasta.py          header-style registry and FASTA/TSV writers
  parse/inputs.py   VCF, ANNOVAR, fusion and SUPPA2 parsers
  build/smallvar.py SNV / InDel / RNA-editing → protein
  build/fusion.py   chimeric transcript construction → protein
  build/splicing.py AS event → isoform protein
scripts/            00 fetch refs, 01 parse, 02 build FASTA, 03 QC
tests/              fixture generator + test suite
results/            manifest, QC tables, FASTA output
logs/               operation logs + provenance JSON
docs/               format specification
```

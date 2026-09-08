# Cross-tool comparison

**Status: pypgatk measured. ProteoDisco not attempted.**

## Result: v2p 260/260, pypgatk 259/260

| tool | version | loci with at least one variant entry | recall |
|---|---|---:|---:|
| **v2p** | 1.0.0 | **260 / 260** | **100%** |
| pypgatk | 0.0.24 | 259 / 260 | 99.6% |

pypgatk misses `chr4:152332595`, a five-base FBXW7 deletion whose VCF
anchor base sits outside the coding exon. v2p missed it too until that was
fixed; the fix is in the history and the regression test is in
`tests/test_pipeline.py`.

Two caveats, unchanged, and a reader cannot interpret the number without
them.

**The two tools were not given the same input.** pypgatk cannot run on a
sites-only VCF: it requires VEP annotation, and VEP is what told it which
transcripts to use. v2p resolved transcripts itself from the unannotated
VCF. That is the only comparison available, since the alternative is not
running pypgatk at all, but it is not a like-for-like test of the same
task - and if anything it favours pypgatk, which is handed the answer to
part of the problem.

**Only one of the two can be scored on protein-change correctness.**
pypgatk's headers carry a locus and a transcript and no protein change,
and its sequences contain internal stop codons - translated transcript
products rather than annotated proteoforms. "Did the right protein change
come out" cannot be asked of it. v2p's 99.2% recall and 99.2% precision on
protein change have no counterpart here.

### Earlier correction

An earlier version of this file reported **99.2% for v2p against 99.6% for
pypgatk**, a loss by one variant. That was wrong, and the fault was in the
measurement. The metric read `LOC=` from the shipped FASTA headers, and
cross-class deduplication merges identical sequences keeping one header,
so a locus whose protein duplicates another variant's disappeared. Two
HRAS variants both encoding `p.F156L` collapsed into one entry and the
second was scored as a miss. Scoring now reads `tables/provenance.jsonl`,
and the underlying provenance loss is fixed at source.

## Getting pypgatk to run at all

pypgatk **0.0.24**, installed from PyPI, given the same truth VCF and
GENCODE v44:

```bash
pypgatk_cli.py vcf-to-proteindb \
    --input_fasta ref/gencode.v44.pc_transcripts.fa.gz \
    --vcf bench_work/inputs/truth_variants.vcf \
    --gene_annotations_gtf ref/gencode.v44.annotation.gtf.gz \
    --output_proteindb pypgatk.fasta --ignore_filters
```

Two failures came first, both of which exited 0:

1. **Sites-only VCF** - it indexed the GTF for twenty minutes, then
   wrote no output file. `vcf-to-proteindb` reads the transcript and
   consequence from a VEP annotation field in INFO, and the truth VCF
   had none.
2. **Plain-gzip transcript FASTA** - `ValueError: Gzipped files are not suitable for indexing, please use
   BGZF`. Biopython cannot random-access a plain `.gz`.

With a VEP-annotated VCF (`scripts/10_vep_annotate.py`) and the
transcripts uncompressed, it translated 335 of 336 variants into
9,253 entries.

Both prerequisites are worth recording, because neither is visible from
pypgatk's feature list and the first one fails **silently**: exit code 0,
no output file, no warning that the annotation field was missing. A
pipeline checking only the exit status would log a successful run that
produced nothing - the exact failure mode v2p's release invariants exist
to catch.

The input difference is the more important one for interpreting the
result above. pypgatk requires a VEP-annotated VCF; v2p takes a sites-only
VCF and does its own transcript resolution and consequence calling. A user
choosing between them is choosing between assembling a VEP annotation
first or not.

## Why the comparison is narrow

Only the small-variant category is comparable at all. No other tool accepts
an editing table, so the ADAR set has no comparator, and no database
builder ingests fusion calls. Those are coverage differences, not accuracy
differences, and claiming victory in a category the competition cannot
enter would be dishonest.

| category | comparable against | why |
|---|---|---|
| small variants | ProteoDisco, pypgatk, customProDB, QUILTS | all four build variant protein databases from VCF |
| splice junctions | ProteoDisco, customProDB, QUILTS | pypgatk does not do junctions |
| A-to-I editing | **nobody** | no other tool accepts an editing table |
| fusions | **nobody** | AGFusion models fusions but emits no search database |
| non-canonical ORFs | pypgatk | v2p gained this in M5; pypgatk had it first |

The honest framing: v2p should be **comparable** on small variants and
**unique** on editing and fusions. If ProteoDisco beats it on the shared
subset, that belongs in the README, not in a drawer.

## The shared input

Use `truth/cosmic_variants.tsv`, the 336 ClinVar rows, converted to VCF.
`scripts/09_benchmark.py` already writes exactly that file — run it once
and keep the intermediate:

```bash
python3 scripts/09_benchmark.py --ref ref/ --workdir bench_work
# bench_work/inputs/truth_variants.vcf is the shared input
```

Score every tool on the same asserted rows (263 rows, 260 distinct
loci), so the comparison is like for like.
`--min-review-status any` gives the all-336 figure as a secondary number.

## ProteoDisco (R / Bioconductor)

```r
install.packages("BiocManager")
BiocManager::install("ProteoDisco")          # record the exact version
packageVersion("ProteoDisco")

library(ProteoDisco)
pd <- generateProteoDiscography(
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    genomeSeqs = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
pd <- importGenomicVariants(pd, files = "bench_work/inputs/truth_variants.vcf")
pd <- incorporateGenomicVariants(pd)
exportProteoDiscography(pd, file = "proteodisco.fasta")
```

Note the annotation mismatch: ProteoDisco is driven from a UCSC TxDb while
v2p uses GENCODE. Different transcript sets mean different canonical
choices, so score on locus and protein change rather than transcript id —
the rule `09_benchmark.py` already applies.

## pypgatk (Python)

```bash
pip install pypgatk                          # record the exact version
pypgatk_cli.py --version

pypgatk_cli.py vcf-to-proteindb \
    --config_file config/vcf_config.yaml \
    --vep_annotated_vcf bench_work/inputs/truth_variants.vcf \
    --input_fasta gencode.v44.pc_transcripts.fa \
    --gene_annotations_gtf ref/gencode.v44.annotation.gtf.gz \
    --output_proteindb pypgatk.fasta
```

pypgatk expects a VEP-annotated VCF and an **uncompressed** transcript
FASTA. Annotate first with `scripts/10_vep_annotate.py`, which uses the
Ensembl REST API and needs no local VEP cache:

```bash
python3 scripts/10_vep_annotate.py --vcf bench_work/inputs/truth_variants.vcf     --out-vcf truth_vep.vcf --out-tsv vep_consequences.tsv
gunzip -k gencode.v44.pc_transcripts.fa.gz
```

Then pass `--annotation_field_name CSQ --consequence_str Consequence
--transcript_str Feature --biotype_str BIOTYPE`, matching the CSQ layout
that script emits.

## customProDB and QUILTS

customProDB is Ensembl-89-era and unmaintained; QUILTS is dormant. Include
them only if a reviewer asks. If they are run, say plainly which reference
release each was given, because neither can be pointed at GENCODE v44.

## Results

| tool | version | locus coverage | protein-change recall | precision |
|---|---|---:|---:|---:|
| v2p | 1.0.0 | **100% (260/260)** | 99.2% | 99.2% |
| pypgatk | 0.0.24 | 99.6% (259/260) | not reportable | not reportable |
| ProteoDisco | — | not attempted | not attempted | not attempted |

pypgatk emits no protein change, so those two columns cannot be filled for
it. ProteoDisco needs R, Bioconductor, BSgenome and a TxDb, and was not
attempted.

Only the v2p row is measured. Fill the others in from real runs, record the
exact versions, and report whatever comes out — including a result that
favours a competitor.

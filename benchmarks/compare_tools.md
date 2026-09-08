# Cross-tool comparison

**Status: pypgatk measured. ProteoDisco not attempted.**

## Result: pypgatk edges v2p on locus coverage

| tool | version | loci with at least one variant entry | recall |
|---|---|---:|---:|
| v2p | 1.0.0 | 258 / 260 | 99.2% |
| **pypgatk** | **0.0.24** | **259 / 260** | **99.6%** |

pypgatk found one locus v2p did not. v2p found none pypgatk did not. One
locus neither produced. **That is a loss, by one variant, and it belongs
here rather than in a drawer.**

Two things must be said alongside it, not to explain the result away but
because a reader cannot interpret the number without them.

**The two tools were not given the same input.** pypgatk cannot run on a
sites-only VCF at all; it requires VEP annotation, and VEP is what told it
which transcripts to use. v2p resolved transcripts itself from the
unannotated VCF. So this compares v2p-doing-its-own-annotation against
pypgatk-given-VEP's-annotation. That is the only comparison available,
because the alternative is not running pypgatk at all, but it is not a
like-for-like test of the same task.

**Only one of the two can be scored on protein-change correctness.**
pypgatk's headers carry a locus and a transcript and no protein change,
and its sequences contain internal stop codons - it emits translated
transcript products rather than annotated proteoforms. So "did the right
protein change come out" cannot be asked of it. v2p's 97.3% recall and
98.1% precision on protein change have no pypgatk counterpart, and it
would be dishonest to present the 99.2%/99.6% row as if it were the same
measurement.

Coverage and correctness are different claims. On coverage of these 260
loci, pypgatk is ahead by one. On protein-change correctness, there is no
comparison to make.

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
| v2p | 1.0.0 | 99.2% (258/260) | 97.3% | 98.1% |
| pypgatk | 0.0.24 | **99.6% (259/260)** | not reportable | not reportable |
| ProteoDisco | — | not attempted | not attempted | not attempted |

pypgatk emits no protein change, so those two columns cannot be filled for
it. ProteoDisco needs R, Bioconductor, BSgenome and a TxDb, and was not
attempted.

Only the v2p row is measured. Fill the others in from real runs, record the
exact versions, and report whatever comes out — including a result that
favours a competitor.

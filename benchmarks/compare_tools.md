# Cross-tool comparison

**Status: pypgatk installed and run; it produced no output, for a reason
worth recording. ProteoDisco not attempted.** No competitor number has
been measured, so v2p's own figures imply nothing about the competition.
Publishing this table with invented numbers would be worse than leaving it
empty.

## What happened when pypgatk was actually run

pypgatk **0.0.24**, installed from PyPI, given the same truth VCF and
GENCODE v44:

```bash
pypgatk_cli.py vcf-to-proteindb \
    --input_fasta ref/gencode.v44.pc_transcripts.fa.gz \
    --vcf bench_work/inputs/truth_variants.vcf \
    --gene_annotations_gtf ref/gencode.v44.annotation.gtf.gz \
    --output_proteindb pypgatk.fasta --ignore_filters
```

It spent about twenty minutes building a gffutils SQLite index of the GTF,
then **exited 0 and wrote no output file at all**.

The cause is not a bug, it is a prerequisite: `vcf-to-proteindb` reads the
transcript id and the consequence from a **VEP annotation field in the VCF
INFO column**. The truth VCF is sites-only — its INFO is `GENE=<symbol>`
and nothing else — so no variant carried an annotation pypgatk could use,
and it had nothing to translate.

Two things follow, and both belong in any comparison:

1. **The tools take different inputs.** pypgatk requires a VEP-annotated
   VCF. v2p takes a sites-only VCF and does its own transcript resolution
   and consequence calling. That is a real difference in what a user has
   to assemble before either tool will run, and it is not visible from
   either tool's feature list.
2. **pypgatk fails silently here.** Exit code 0, no file, no warning that
   the annotation field was missing. A pipeline that checked only the exit
   status would record a successful run that produced nothing. This is
   exactly the failure mode v2p's release invariants exist to catch, and
   it is worth stating plainly rather than treating as an implementation
   detail.

To finish this comparison, the truth VCF must be VEP-annotated first —
VEP plus its GRCh38 cache is roughly 25 GB, which is why it has not been
done here. Once annotated, rerun the command above and fill in the table.

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

Score every tool on the same 263 asserted rows, with the same
locus-plus-protein-change matching, so the comparison is like for like.
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

pypgatk expects a VEP-annotated VCF. The truth VCF is not annotated, so
either run VEP over it first or use pypgatk's own Ensembl route. Whichever
is chosen must be recorded here, because it changes what is being measured.

## customProDB and QUILTS

customProDB is Ensembl-89-era and unmaintained; QUILTS is dormant. Include
them only if a reviewer asks. If they are run, say plainly which reference
release each was given, because neither can be pointed at GENCODE v44.

## Results

| tool | version | small-variant recall | precision | notes |
|---|---|---|---|---|
| v2p | 1.0.0 | 97.3% | 98.1% | 263 asserted rows, `09_benchmark.py` |
| pypgatk | 0.0.24 | not measured | not measured | ran, exited 0, produced no output: needs a VEP-annotated VCF |
| ProteoDisco | — | not attempted | not attempted | needs R, Bioconductor, BSgenome + TxDb |

Only the v2p row is measured. Fill the others in from real runs, record the
exact versions, and report whatever comes out — including a result that
favours a competitor.

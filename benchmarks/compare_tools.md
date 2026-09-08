# Cross-tool comparison

**Status: specified, not run.** Everything below is procedure. No
competitor has been installed or executed, so this document contains no
competitor results, and v2p's own figures imply nothing about the
competition. Filling the results table in is the remaining work;
publishing it with invented numbers would be worse than leaving it empty.

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
| ProteoDisco | — | not run | not run | |
| pypgatk | — | not run | not run | |

Only the v2p row is measured. Fill the others in from real runs, record the
exact versions, and report whatever comes out.

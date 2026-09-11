# Benchmark truth sets

Truth data, and the runner that scores against it.

```bash
python3 scripts/09_benchmark.py --ref ref/                          # 263 rows
python3 scripts/09_benchmark.py --ref ref/ --min-review-status any  # all 336
```

Latest result: **99.2% recall, 99.2% precision** on the 263 small-variant
rows with stated assertion criteria; 99.4% and 99.4% across all 336 -
slightly higher, so excluding the unasserted rows is not flattering the
figure.
A-to-I editing 5 of 5. Full report, including every disagreement, in
`results/benchmark.md`.

This scores v2p against a truth set. It is **not** a comparison against
ProteoDisco or pypgatk — neither has been run — so it supports no claim to
beat another tool.

Two things worth knowing about how it scores. Matching is on locus plus
protein change, not transcript id, because the truth set records RefSeq
transcripts while the pipeline works in GENCODE, and mapping between them
would add a second source of error to the measurement. And a frameshift
written `p.M862fs` by ClinVar is accepted as matching `p.M862Ifs*4` from
v2p — the same event at different precision. A different residue,
position, or consequence class is never forgiven.

## `truth/cosmic_variants.tsv` — 336 small variants

Retrieved from NCBI ClinVar via E-utilities by
`build_truth_from_clinvar.py`. GRCh38 throughout.

Established per row: the ClinVar VCV accession and its version, the
SPDI-derived coordinate, and the REF allele, which was re-read from
`ref/GRCh38.primary_assembly.genome.fa` and confirmed to match — 336 of
336, checked independently of the builder that wrote the file.

**Not** established: that each expected protein consequence is right in the
sense a benchmark needs. ClinVar's own annotation is taken at face value.
The `clinvar_review_status` column is the thing to filter on — 8 rows are
"reviewed by expert panel", 26 "criteria provided, multiple submitters, no
conflicts", 229 "criteria provided, single submitter", and **73 carry "no
assertion criteria provided"**. Treat that last group as unreviewed.

The filename says `cosmic` for historical continuity; the contents are
ClinVar. COSMIC requires a licence for bulk
download.

## `truth/editing_sites.tsv` — 5 ADAR recoding sites

The canonical A-to-I recoding events: GRIA2 Q607R, NEIL1 K242R, BLCAP Y2C,
CDK13 Q103R, COG3 I635V. No competing tool accepts an editing table, so
this category cannot be run against customProDB, QUILTS, ProteoDisco or
pypgatk at all.

Coordinates here are **derived, not quoted**. For each published protein
change the codon is located in the GENCODE v44 representative transcript,
and the row is written only if a single A→G in that codon reproduces the
published substitution and the genome base matches the expected strand.
The `verification_method` column records the codon and the substitution
for each row, so the derivation can be rechecked by hand.

This matters: four of these five coordinates were first written from
memory and were **wrong**. The derivation caught all four. Do not
reintroduce hand-typed coordinates here.

## Checks the README does not have room for

The README's "Measured" table is capped at six rows. The rest of what is
verified:

| check | result | where |
|---|---|---|
| ClinVar, all 336 rows including unasserted | 99.4% recall and precision | `results/benchmark.md` |
| A-to-I editing positive controls | 5/5 | `truth/editing_sites.tsv`, `results/benchmark.md` |
| GENCODE translation agreement | 100% on a single-transcript build | the `--min-agreement` gate, and `qc/` in any release |
| combinatorial entries kept vs dropped | 9 kept, 40 dropped for adding no peptide | [USAGE](../docs/USAGE.md#combinatorial-proteoforms) |

## ProVar comparison

`provar_comparison.md` holds the table; `run_provar.sh` rebuilds it from
scratch. ProVar is the fair comparator — same input shape, same
sample-level scope, same authors as ProHap — and it was run on the same
HCC1395 VCF against Ensembl 110, the release GENCODE v44 is built from.

### What ProVar represents and v2p does not

Eight variants, and every one was investigated rather than waved away:

| variant | why v2p produced nothing | verdict |
|---|---|---|
| `19:37692270G>A` | `below_min_length`: the product is under 8 aa. ProVar's own row for it reads `-160:Q>*` — a stop upstream of the annotated start | policy difference; a 7-residue entry is not searchable |
| 7 deletions, 28–46 bp | each begins inside a coding exon and runs past its boundary, removing a splice site | deliberate refusal |

The deletions are the interesting case. ProVar translates the exonic
part and records the splice damage in a `splice_site_affected` column
(values of 2, 4 and 17 bases for these). v2p declines, because trimming
a deletion to the exon boundary asserts a spliced product nobody has
observed — the 37 bp deletion at `1:155086199` loses only two residues
that way, which tells you most of it was never exonic.

That refusal stands. What did not stand is how it was recorded: all
seven were filed as `not_in_coding_exon`, which is wrong twice over —
they *are* in coding exons, and the reason they were skipped, that they
destroy a splice site, is the part a reader wants. They now get their
own disposition outcome, `crosses_splice_junction`, with the span and
the reason. Five HCC1395 variants move into it, the database is
unchanged, and `test_pipeline.py` carries the regression test.

## The search for a cross-evidence dataset

v2p's one unique claim is combining an RNA editing site with a DNA
variant on one transcript. Demonstrating it needs a sample whose calls
contain both, co-located. HCC1395 does not, and that was measured rather
than assumed:

| check | count |
|---|---|
| transcripts carrying both an editing site and a DNA variant | 40 |
| of those, transcripts where **both** recode | **0** |
| genes with a recoding edit *and* a recoding DNA variant | **0** |

29 genes carry a recoding edit, 278 carry a recoding DNA variant, and
the sets do not intersect. The 40 transcripts that carry both have the
editing site in a UTR or at a synonymous position, so no combined
proteoform arises.

Public alternatives were searched for. SEQC2 publishes WGS, WES and
RNA-seq for HCC1395 but no consensus A-to-I call set. REDIportal
publishes 4,388 nonsynonymous recoding sites, but its bulk download is
behind a form and the direct paths tried all returned 404. Nothing found
offers both call types, on one sample, in a form that can be downloaded
and re-run.

So `cross_evidence_case.md` is **constructed**, and says so: a real
recoding editing site from HCC1395 plus a synthetic SNV placed in the
same tryptic peptide. It shows the mechanism, not how often such
proteoforms occur.

## Not present

`truth/fusions.tsv` is not included: the published breakpoints for BCR–ABL1, EML4–ALK and
TMPRSS2–ERG could not be verified to the standard of the two files above,
and an unverified fusion truth set is worse than none. The three fusions
in `examples/fusions.csv` are illustrative input, not truth data.

A ProteoDisco comparison is not included: its R/Bioconductor environment
could not be installed here, and `run_proteodisco.R` is the driver for
anyone who can. The pypgatk comparison is in `compare_tools.md`.

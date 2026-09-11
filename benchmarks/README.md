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

## Not present

`truth/fusions.tsv` is not included: the published breakpoints for BCR–ABL1, EML4–ALK and
TMPRSS2–ERG could not be verified to the standard of the two files above,
and an unverified fusion truth set is worse than none. The three fusions
in `examples/fusions.csv` are illustrative input, not truth data.

A ProteoDisco comparison is not included: its R/Bioconductor environment
could not be installed here, and `run_proteodisco.R` is the driver for
anyone who can. The pypgatk comparison is in `compare_tools.md`.

# Benchmark truth sets

Truth data only. **There is no benchmark runner yet**, so no precision or
recall figure has been computed from these files and v2p makes no accuracy
claim. Read that before citing anything here.

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

The filename says `cosmic` for continuity with the `docs/BUILD_SPEC.md`
contract; the contents are ClinVar. COSMIC requires a licence for bulk
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

`truth/fusions.tsv` is specified in `docs/BUILD_SPEC.md` but is not
included: the published breakpoints for BCR–ABL1, EML4–ALK and
TMPRSS2–ERG could not be verified to the standard of the two files above,
and an unverified fusion truth set is worse than none. The three fusions
in `examples/fusions.csv` are illustrative input, not truth data.

`scripts/09_benchmark.py` — the runner that would consume these files and
report per-category precision and recall — is also not written.

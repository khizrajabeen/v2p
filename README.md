# v2p — variant calls to protein sequences

v2p turns genomic and transcriptomic variant calls into an amino-acid
sequence database for mass-spectrometry search, and records what happened
to every input variant. It reads somatic SNVs and indels, A-to-I RNA
editing tables, gene fusion calls and alternative-splicing events, and
emits one coherent FASTA with a shared header vocabulary. It exists
because no other tool covers those four evidence types together: RNA
editing produces proteoforms with no genomic basis, so they are invisible
to any DNA-driven pipeline, and fusions and splicing normally need
separate pipelines whose outputs nobody reconciles. Every variant that
does *not* produce a protein gets a recorded reason, because silence is
where wrong answers hide.

## Where it sits

| tool | year | language | small variants | splice junctions | RNA editing | fusions | non-canonical ORFs | status |
|---|---|---|---|---|---|---|---|---|
| customProDB | 2013 | R | yes | yes | no | no | no | unmaintained |
| QUILTS | 2016 | Python | yes | yes | no | no | partial | dormant |
| ProteoDisco | 2021 | R/Bioconductor | yes | yes | no | no | no | maintained |
| pypgatk / pgdb | 2021 | Python | yes | no | no | no | **yes** (3-frame) | maintained |
| **v2p** | — | Python | yes | yes | **yes** | **yes** | no | this repo |

The claim worth making, and only this one: the first tool to build a single
protein database from DNA variants, RNA editing, fusions and splicing
together, with every input variant's fate recorded. Not "more accurate" —
that needs a published benchmark this repo does not yet have. Not "better
than VEP" — VEP is not a database builder. Three-frame translation of
non-coding transcripts is a real gap where pypgatk is ahead.

## Install

```bash
pip install .
```

Python 3.10–3.13. The only runtime dependency is `pyfaidx`.

Then fetch the references once (~18 GB, GENCODE v44 + GRCh38):

```bash
bash scripts/00_fetch_references.sh ref/
```

## Quickstart

`examples/` holds a small synthetic input folder that ships with the repo,
so this works immediately after cloning:

```bash
v2p detect examples/
```

```
gene fusion calls                        fusions.csv
RNA editing sites                        rna_editing.txt
somatic SNV / InDel calls                small_variants.vcf
alternative splicing events              splicing.csv
genome build                             GRCh38
```

With the references in place, the whole conversion is one command:

```bash
v2p run examples/ --ref ref/ --outdir out/ --name EXAMPLE --logdir logs/
```

It prints the invariant report and exits non-zero if the release
contradicts itself. Keep `--logdir` outside `--outdir`: `MANIFEST.txt` is
written last, so anything added to the release afterwards makes it stale.

The example calls are synthetic calls at real coordinates — nothing is
redistributed from the HCC1395 package. `examples/make_examples.py`
regenerates them: substitutions are placed in the CDS of named genes with
the reference allele read from the genome, and the ADAR sites are derived
from the annotation rather than quoted (see Validation).

## Validation

- **GENCODE agreement 100%** on the single-transcript build — reference
  proteins translated by this pipeline compared against GENCODE's own
  translations of the same transcript ids. 97.3% on the all-transcript
  build; the difference is `cds_end_NF` transcripts, where GENCODE
  truncates at the annotated CDS end and v2p reads to the first stop.
  This is the correctness gate, not UniProt: UniProt and GENCODE choose
  canonical isoforms independently and often disagree on the start codon.
- **Nonsynonymous:synonymous ratio 2.75**, within the expected range.
- **ADAR positive controls recovered**: GRIA2 Q607R, NEIL1 K242R,
  BLCAP Y2C, CDK13 Q103R, COG3 I635V. Their coordinates in
  `benchmarks/truth/editing_sites.tsv` are *derived*, not quoted: for each
  published protein change the codon is located in the GENCODE v44
  representative transcript, and the row is written only if a single A→G
  in that codon reproduces the published substitution and the genome base
  matches the expected strand.
- **200 tests**, all offline, no reference download, seconds to run:
  93 pipeline, 80 release-invariant, 27 config.
- **Nine release invariants** run before `v2p run` reports success, and
  any error-severity violation exits non-zero. On the HCC1395 dataset the
  release reports zero errors.
- **Reproducibility**: two runs from one config produce a byte-identical
  release tree, verified with `diff -r` over the entire output including
  `MANIFEST.txt`.

```bash
make test              # 200 assertions, no reference needed
make reproducibility   # double-run byte-identity, needs the references
```

## How it works

Variants are applied to the **mature transcript**, not to genomic
sequence. That is what lets one code path serve DNA variants and RNA
editing: an A-to-I edit is indistinguishable from a genomic A>G at
transcript level. Translation runs from the annotated start codon to the
first stop; frameshifts translate past the annotated stop into the 3′ UTR.

Genomic coordinates are 1-based inclusive; transcript and CDS offsets are
0-based from the 5′ end in transcription direction. Every conversion goes
through `annotation.py` so the convention lives in one place. The genetic
code is hard-coded in `seqops.py` on purpose — a library update must not
be able to change a translation silently.

Transcript selection ties break deterministically: MANE_Select →
Ensembl_canonical → basic → longest CDS → longest transcript →
lexicographic id.

## Limitations and what is unverified

**Not implemented.** Non-canonical ORFs (three-frame translation of
lncRNAs and pseudogenes) — the one capability where pypgatk is genuinely
ahead. Species independence: `_HUMAN` entry names, `OS=Homo sapiens
OX=9606` and GRCh38 contig lengths are still hard-coded, so this is a
human-only tool today. There is no benchmark runner: `benchmarks/truth/`
holds the truth data but no script consumes it, so no precision or recall
figure has been computed and no accuracy claim is made.

**Unverified.** The GitHub Actions workflow has never executed — it will
run on first push and prove itself or not. The `Dockerfile` has never been
built; treat it as a starting point, not a tested artefact. There is no
conda recipe: one was written and then deleted unbuilt, rather than
shipped untested.

**Truth set.** `benchmarks/truth/cosmic_variants.tsv` holds 336 ClinVar
variants, each with its VCV accession and version, GRCh38 throughout, and
every REF allele independently re-read from the genome and confirmed. What
is *not* established is that each expected protein consequence is right in
the sense a benchmark needs — ClinVar's own annotation is taken at face
value, and 73 of the 336 rows carry "no assertion criteria provided". Use
it accordingly.

**Bugs found and fixed during development.** Every one produced
plausible-looking *wrong output* rather than an error, and none was caught
by looking at the sequences. Each has a regression test.

| bug | effect |
|---|---|
| transcripts looked up by gene symbol only | 33,221 sites-only VCF variants produced zero proteins |
| minus-strand AF/AL events dropped | 234 of 1,027 splicing events silently skipped |
| GTF frame column parsed, never used | 8,260 `cds_start_NF` transcripts translated out of frame |
| selenocysteine annotation ignored | all 25 selenoproteins truncated at their first UGA |
| mitochondrial code not applied | chrM would read through its real stops |
| UTR variants classified as synonymous | inverted the nonsyn:syn ratio; 16 indels shipped as false frameshifts |
| ambiguous gene symbols resolved by name | silent false negatives on paralogous loci |
| double-counted wild-types across files | per-type files did not sum to the combined file |
| `v2p run --ref` ignored the reference directory | the documented invocation always refused without `--force` |
| `--logdir` reached only the invariant logger | one run's provenance scattered across two directories |

An eleventh, found while writing the examples: four of the five ADAR
coordinates first written from memory were wrong. They are now derived
from the annotation, and the derivation fails loudly rather than guessing.

## Layout

```
src/v2p/        cli.py config.py discover.py invariants.py
                annotation.py seqops.py nmd.py peptides.py validate.py fasta.py
                provenance.py  parse/  build/  stages/
tests/          test_pipeline.py test_invariants.py test_config.py
docs/           BUILD_SPEC.md TOOL_DESIGN.md FORMAT_SPEC.md
examples/       synthetic input, works after clone
benchmarks/     truth data (no runner yet)
config/         params.yaml — generated, and tested to load
```

The numbered pipeline stages live in `src/v2p/stages/` so a pip-installed
v2p can find them. Each is a separate process so it records its own
provenance JSON — inputs with checksums, parameters, counters, environment.

## Citation

A paper is not yet written. If you use v2p before then, please cite the
repository and the exact commit:

```
v2p: variant calls to protein sequences. https://github.com/<owner>/v2p
Version 1.0.0, commit <sha>.
```

See `CITATION.cff`.

## Licence

**Not yet chosen.** The copyright holder is being confirmed; until a
`LICENSE` file lands, no licence is granted and the default of "all rights
reserved" applies. `pyproject.toml` and `CITATION.cff` carry TODOs marking
the same gap.

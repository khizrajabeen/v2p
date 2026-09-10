# v2p — variant calls to protein sequences

[![tests](https://github.com/khizrajabeen/v2p/actions/workflows/tests.yml/badge.svg)](https://github.com/khizrajabeen/v2p/actions/workflows/tests.yml)
[![reference audit](https://github.com/khizrajabeen/v2p/actions/workflows/reference-audit.yml/badge.svg)](https://github.com/khizrajabeen/v2p/actions/workflows/reference-audit.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![assemblies](https://img.shields.io/badge/assemblies-GRCh38%20%7C%20GRCh37%20%7C%20T2T-orange)](#installation)

**v2p** builds mass-spectrometry search databases from variant calls, and
records what happened to **every** input variant — including the ones that
produce no protein. It reads DNA variants, RNA editing, gene fusions and
alternative splicing into one FASTA with a shared header vocabulary.

## Key features

- **Four evidence types in one database** — somatic SNVs/indels, A-to-I
  RNA editing, gene fusions, alternative splicing. No other tool covers
  this combination.
- **RNA editing as a first-class input** — ADAR recoding produces
  proteoforms with no genomic basis, invisible to any DNA-driven pipeline.
- **Per-variant disposition** — every input gets an outcome: protein
  built, synonymous, UTR, intronic, no coding transcript, error.
- **Nine release invariants** — mechanical checks that run before success
  is reported; any error-severity violation exits non-zero.
- **Provenance graph** — one record per output sequence naming every
  variant that produced it, and one row per novel peptide naming the
  variants that explain it.
- **Content-based input detection** — point it at a folder; no filename
  convention required.
- **Reference-interface auditing** — verifies the GTF, genome and proteome
  are being *read* correctly, separately from whether output looks right.
- **Byte-identical reproducibility** — two runs from one config produce
  the same release, verified over the whole tree.
- **Combinatorial proteoforms** — when two or more variants land on one
  transcript, emit the protein carrying *all* of them. A tryptic peptide
  spanning two variants is in neither single-variant entry nor the
  reference, so a database built one variant at a time cannot identify it
  at any FDR; [ProHap](https://doi.org/10.1038/s41592-024-02506-0)
  measured this at 12.4% of substitutions for common germline haplotypes.
  Two variants can even share a codon: in HCC1395, FKTN carries
  `p.[D225K]`, a lysine no single-variant entry produces, which also
  creates a tryptic cleavage site. Each entry must contribute a peptide no
  single-variant entry has, or it is dropped rather than inflating the
  database. Phase is honoured where the caller reports `GT`/`PS` —
  variants on opposite haplotypes are never combined — and entries are
  marked `unphased` where it is absent. Opt-in with `--combine-variants`;
  [how v2p compares](#how-v2p-compares) says what is and is not new here.
- **Non-canonical ORFs** — three-frame translation of lncRNAs,
  pseudogenes, and 5′/3′ UTRs of coding transcripts (uORFs and dORFs),
  with the reading frame recorded per entry and a configurable length
  floor. Entries identical to a reference protein are dropped rather than
  duplicating the database. Opt-in with `--include-noncanonical`, since
  three-frame translation multiplies database size.
- **Multi-assembly** — GRCh38, GRCh37 and T2T-CHM13 all fetchable and
  supported, with the build asserted by contig length before translation
  starts.
- **Multi-species** — human, mouse, rat and zebrafish ship; others are one
  YAML file.
- **Search-engine compatible** — four header conventions (UniProt, PEFF,
  pVAC, descriptive) and three decoy strategies.

## How v2p compares

Combining co-occurring variants is not new, and v2p did not invent it.
It is the central purpose of
[ProHap](https://doi.org/10.1038/s41592-024-02506-0) (*Nature Methods*
2024), which builds protein haplotypes "using observed combinations of
alleles in each transcript" from **phased genotype data**. Its companion
**ProVar** takes sample-level VCFs and "considers each allele
independently" — one sequence per variant, no combination. v2p works in
ProVar's scope, one sample, and combines there.

| | ProHap | ProVar | v2p |
|---|---|---|---|
| scope | phased genotype panels | one sample | one sample |
| combines co-occurring variants | **yes**, phased germline | no | yes |
| RNA editing | no | no | **yes** |
| fusions | no | no | **yes** |
| splicing | transcript-level | no | **yes** |
| install | Snakemake + Conda, ~1 TB for full 1000G | Snakemake + Conda | pip, one dependency |
| output audit trail | peptide annotator | — | disposition + provenance graph |

**The one claim that is v2p's alone:** it combines co-occurring variants
*across evidence types*. A somatic SNV co-occurring with an A-to-I edit on
one transcript cannot be represented by any haplotype panel — RNA editing
is not in the genome and never appears in a VCF of genotypes. ProHap
cannot see that proteoform in principle, not by omission.

## Installation

```bash
pip install .
```

Python 3.10–3.13. The only runtime dependency is `pyfaidx`.

Then fetch a reference once:

```bash
bash scripts/00_fetch_references.sh ref/                     # GRCh38, the default
bash scripts/00_fetch_references.sh ref37/  --assembly GRCh37
bash scripts/00_fetch_references.sh reft2t/ --assembly T2T
```

Each downloads the genome, annotation and — where the source publishes
them — the annotation's own protein translations, then asserts the build
by contig length before you spend an hour translating against the wrong
one. `GENCODE_RELEASE=45 bash scripts/00_fetch_references.sh ref/` picks a
different release.

## Quick start

`examples/` ships with the repo, so this works immediately after cloning:

```bash
# 1. see what is in a folder — writes nothing
v2p detect examples/ --ref ref/

# 2. convert
v2p run examples/ --ref ref/ --outdir out/ --name MYSAMPLE --logdir logs/

# 3. check the reference interface
v2p audit --ref ref/
```

## Commands

### Main interface

| command | description |
|---|---|
| `v2p detect <folder>` | identify inputs by content and print the run plan |
| `v2p run <folder>` | the whole conversion, ending in the invariant report |
| `v2p audit --ref <dir>` | verify the GTF, genome and proteome are read correctly |

### Choosing a reference

Any reference release or assembly works — nothing is pinned in code.
Either point at a directory and let detection sort it out, or name each
file:

| option | description |
|---|---|
| `--ref <dir>` | directory holding genome, annotation and proteome |
| `--genome <fa>` | reference genome FASTA, explicit |
| `--annotation <gtf>` | GTF/GFF annotation, explicit |
| `--proteome <fa>` | reference proteome FASTA, explicit |
| `--translations <fa>` | the annotation source's own protein translations |
| `--species <name\|yaml>` | `human`, `human_grch37`, `human_t2t`, `mouse`, `rat`, `zebrafish`, or a YAML path |

```bash
# a different assembly or release, no --ref at all
v2p run calls/ \
    --genome       /refs/GRCh38.primary_assembly.genome.fa \
    --annotation   /refs/gencode.v44.annotation.gtf.gz \
    --proteome     /refs/uniprot_human_SP.fasta \
    --translations /refs/gencode.v44.pc_translations.fa.gz \
    --outdir out/
```

Supply `--translations` whenever you can. It enables the GENCODE
translation check, the gate that actually decides whether the conversion
is right. Without it, validation falls back to comparing against UniProt
canonical, which disagrees on start codons often enough to fail the
threshold for reasons that are not errors.

### Naming inputs explicitly

Detection can be bypassed per evidence type:

| option | description |
|---|---|
| `--small-variants <vcf>` | somatic SNV/InDel calls |
| `--rna-editing <txt>` | A-to-I editing table (ANNOVAR-style) |
| `--fusion-calls <csv>` | gene fusion calls |
| `--splicing <csv>` | alternative splicing events (SUPPA2 ids) |

### Output control

| option | description |
|---|---|
| `--header-style uniprot\|peff\|pvac\|descriptive` | re-emit without re-translating |
| `--decoys none\|reverse\|pseudo_reverse\|shuffle` | target-decoy strategy |
| `--transcript-mode all\|representative` | every transcript, or one per gene |
| `--drop-unchanged` | drop synonymous and UTR variants (kept by default) |
| `--combine-variants` | emit proteins carrying all co-occurring variants on a haplotype |
| `--include-noncanonical` | three-frame translate non-coding transcripts |
| `--split-by-type` / `--no-split` | per-variant-type FASTA files |

### Reproducibility

| option | description |
|---|---|
| `--config <yaml>` | run from a config file |
| `--write-config <yaml>` | emit the config a command line implies |
| `--min-agreement <f>` | translation-agreement gate, default 0.90 |
| `--skip-invariant <ID>` | turn off one release check, recorded in provenance |

### Analysis scripts

| script | description |
|---|---|
| `scripts/09_benchmark.py` | score against the truth sets, per category |
| `scripts/10_vep_annotate.py` | annotate a VCF via Ensembl VEP REST (no local cache) |
| `benchmarks/run_proteodisco.R` | run ProteoDisco on the same reference |

## Supported inputs

| type | recognised by | required fields |
|---|---|---|
| small variants | `##fileformat=VCF` header | CHROM, POS, REF, ALT |
| RNA editing | tab-separated with `Func.refGene`/`Func.ensGene` | Chr, Start, Ref, Alt |
| gene fusion | two gene columns and two `chr:pos` columns | gene1, gene2, breakpoint1, breakpoint2 |
| alternative splicing | SUPPA2 event ids | event id, optionally dPSI |
| confidence regions | BED, three or more columns | chrom, start, end |

## Output

| file | what it is |
|---|---|
| `<name>.target.fasta` | the database — search this |
| `<name>.target_decoy.fasta` | the same, with decoys appended |
| `<name>.entries.tsv` | one row per sequence: class, gene, transcript, protein change, peptide counts |
| `<name>.summary_by_class.tsv` | composition per variant type |
| `by_class/` | one FASTA per variant type, so each can be searched separately |
| `tables/provenance.jsonl` | one record per sequence, naming every variant that produced it |
| `tables/peptide_provenance.tsv` | one row per novel peptide, and whether it spans the variant residue |
| `_work/tables/disposition.tsv` | one row per **input** variant and its fate |
| `qc/` | validation and per-class recovery reports |
| `MANIFEST.txt` | SHA-256 of every file; invariant I8 verifies it against disk |

`disposition.tsv` is the one people miss. Every input variant appears in
it exactly once, including the ones that produced nothing, with the
reason — intronic, no coding transcript, synonymous. It is the only way
to state a recovery rate honestly.

## Measured

| check | result |
|---|---|
| ClinVar benchmark, 263 asserted rows | **99.2% recall, 99.2% precision** |
| same, all 336 rows including unasserted | 99.4% / 99.4% |
| Ensembl VEP consequence agreement | **334/334 (100%)** on shared transcripts |
| A-to-I editing positive controls | **5/5** |
| GENCODE translation agreement | **100%** (single-transcript build) |
| vs pypgatk 0.0.24, locus coverage | **260/260** against 259/260 |
| combinatorial entries on HCC1395 | **9 kept**, 40 dropped as adding no peptide |
| tests | **330**, offline, seconds |

```bash
make test                                   # 330 assertions, no reference needed
make reproducibility                        # two runs from one config, byte-identical
python3 scripts/09_benchmark.py --ref ref/  # the benchmark
```

Details, including every disagreement and the caveats on the pypgatk
comparison, in [benchmarks/](benchmarks/).

## How it works

Variants apply to the **mature transcript**, not to genomic sequence.
That is what lets one code path serve DNA variants and RNA editing: an
A-to-I edit is indistinguishable from a genomic A>G at transcript level.
Translation runs from the annotated start to the first stop; frameshifts
read into the 3′ UTR.

Genomic coordinates are 1-based inclusive; transcript and CDS offsets are
0-based from the 5′ end. The genetic code is hard-coded in `seqops.py` so
a library update cannot change a translation silently. Transcript ties
break deterministically: MANE_Select → Ensembl_canonical → basic →
longest CDS → longest transcript → id.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — running your own data
- [docs/FORMAT_SPEC.md](docs/FORMAT_SPEC.md) — header grammar, all four styles

## Licence

MIT — see [LICENSE](LICENSE).

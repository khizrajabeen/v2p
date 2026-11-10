# v2p — variant calls to protein sequences

[![tests](https://github.com/khizrajabeen/v2p/actions/workflows/tests.yml/badge.svg)](https://github.com/khizrajabeen/v2p/actions/workflows/tests.yml)
[![reference audit](https://github.com/khizrajabeen/v2p/actions/workflows/reference-audit.yml/badge.svg)](https://github.com/khizrajabeen/v2p/actions/workflows/reference-audit.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![assemblies](https://img.shields.io/badge/assemblies-GRCh38%20%7C%20GRCh37%20%7C%20T2T-orange)](docs/USAGE.md#step-1--put-the-reference-in-place-once)

Builds mass-spectrometry search databases from variant calls. Reads DNA
variants, RNA editing, gene fusions and alternative splicing into one
FASTA, and records an outcome for every input variant — including the
ones that produce no protein.

## Install

```bash
pip install .                                  # Python 3.10-3.13, one dependency
bash scripts/00_fetch_references.sh ref/       # GRCh38 + GENCODE, once
```

## Use

```bash
v2p detect examples/ --ref ref/                # what is in a folder; writes nothing
v2p run    examples/ --ref ref/ --outdir out/ --name DEMO --logdir logs/
```

Inputs are identified by content, so no naming convention is required. A
sites-only VCF on its own is a valid input. Every option, output file and
failure message is in [docs/USAGE.md](docs/USAGE.md).

## Features

- **Four evidence types in one database** — somatic SNVs and indels,
  A-to-I editing, gene fusions, alternative splicing.
- **Per-variant disposition** — every input gets a recorded outcome, so a
  recovery rate can be stated rather than assumed.
- **[Combinatorial proteoforms](docs/USAGE.md#combinatorial-proteoforms)** —
  the protein carrying all co-occurring variants, including across
  evidence types, and **phase honoured where the caller reports it**:
  variants the caller placed on opposite haplotypes are never combined,
  and combinations without phase are marked as the hypotheses they are.
- **[Non-canonical ORFs](docs/USAGE.md#about---include-noncanonical)** —
  three-frame translation of lncRNAs, pseudogenes and UTRs. Opt-in.
- **Nine release invariants** — mechanical self-checks before success is
  reported; any error-severity violation exits non-zero.
- **Provenance** — one record per sequence naming every variant behind it,
  one row per novel peptide.
- **Byte-identical reproducibility** — two runs from one config produce
  the same release, verified over the whole tree.
- **Multi-assembly and multi-species** — GRCh38, GRCh37, T2T-CHM13;
  human, mouse, rat, zebrafish.

## Measured

| check | result |
|---|---|
| ClinVar benchmark, 263 asserted rows | 99.2% recall and precision |
| Ensembl VEP consequence agreement | 334/334 on shared transcripts |
| GENCODE translation agreement | 99.4% exact on 4,041 transcripts |
| vs pypgatk 0.0.24, locus coverage | 260/260 against 259/260 |
| vs ProVar, variants represented | 2,865 against 911 |
| tests | 369, offline, seconds |

```bash
make test && make reproducibility
```

Every disagreement and how to rerun each comparison:
[benchmarks/](benchmarks/).

## Scope and limits

Neither combining co-occurring variants nor combining across evidence
types is new. [ProHap](https://doi.org/10.1038/s41592-024-02506-0)
(*Nature Methods* 2024) combines phased germline haplotypes;
[moPepGen](https://doi.org/10.1038/s41587-025-02701-0) (*Nature
Biotechnology* 2025) covers SNVs, indels, fusions, splicing, RNA editing
and circRNA, and states it generates peptides for "any combination" of
them. v2p overlaps that scope and does not claim to have got there first.

What v2p does that those papers do not describe:

- **Phase is a constraint, not a footnote.** Two heterozygous variants on
  opposite alleles sit on no single molecule, so that proteoform does not
  exist. v2p reads `GT`/`PS`, refuses those combinations, and marks the
  rest `unphased` so a hypothesis is never mistaken for an observation.
  moPepGen's paper does not discuss phase or haplotypes.
- **Every input variant gets a recorded outcome**, including the ones
  producing nothing, so a recovery rate is measured rather than assumed.
- **The release checks itself** — nine invariants, exiting non-zero on
  any error-severity violation.

Remaining limits:

- Cross-evidence combination is proven by test and by a
  [constructed case](benchmarks/cross_evidence_case.md), not yet by an
  observation in real data.
- v2p has not been benchmarked head to head against moPepGen.
- Where a caller reports no phase, co-occurrence is assumed rather than
  known. Entries say which in `PHASE=`; `--no-allow-unphased` drops them.
- Only genetic-code tables 1 and 2 are implemented. A species needing
  another is refused, not silently mistranslated.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — installation, options, output, troubleshooting
- [docs/FORMAT_SPEC.md](docs/FORMAT_SPEC.md) — header grammar, all four styles
- [benchmarks/](benchmarks/) — truth sets, tool comparisons
- [LICENSE](LICENSE) — MIT

# v2p — variant calls to protein sequences

[![tests](https://github.com/khizrajabeen/v2p/actions/workflows/tests.yml/badge.svg)](https://github.com/khizrajabeen/v2p/actions/workflows/tests.yml)
[![reference audit](https://github.com/khizrajabeen/v2p/actions/workflows/reference-audit.yml/badge.svg)](https://github.com/khizrajabeen/v2p/actions/workflows/reference-audit.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![assemblies](https://img.shields.io/badge/assemblies-GRCh38%20%7C%20GRCh37%20%7C%20T2T-orange)](docs/USAGE.md#step-1--put-the-reference-in-place-once)

**v2p** builds mass-spectrometry search databases from variant calls,
reading DNA variants, RNA editing, gene fusions and alternative splicing
into one FASTA. It records what happened to **every** input variant,
including the ones that produce no protein.

```bash
pip install .                                  # Python 3.10-3.13, one dependency
bash scripts/00_fetch_references.sh ref/       # GRCh38 + GENCODE, once
v2p detect examples/ --ref ref/                # what is in a folder; writes nothing
v2p run    examples/ --ref ref/ --outdir out/ --name DEMO --logdir logs/
```

Every option, output file and failure message is in
[docs/USAGE.md](docs/USAGE.md).

## Why another one

Combining co-occurring variants is not new, and v2p did not invent it.
[ProHap](https://doi.org/10.1038/s41592-024-02506-0) (*Nature Methods*
2024) builds protein haplotypes "using observed combinations of alleles in
each transcript" from phased genotype panels; its companion **ProVar**
"considers each allele independently". v2p works in ProVar's scope — one
sample — and combines there.

| | ProHap | ProVar | pypgatk | v2p |
|---|---|---|---|---|
| scope | phased genotype panels | one sample | one sample | one sample |
| input | phased VCF | VCF | VEP-annotated VCF | a sites-only VCF is enough |
| combines co-occurring variants | **yes**, phased germline | no | no | yes |
| RNA editing | no | no | no | **yes** |
| fusions | no | no | no | **yes** |
| splicing | transcript-level | no | no | **yes** |
| non-canonical ORFs | no | no | yes | yes |
| install | Snakemake + Conda, ~1 TB for full 1000G | Snakemake + Conda | pip | pip |
| per-variant audit trail | peptide annotator | — | — | disposition + provenance graph |

**The one claim that is v2p's alone:** it combines co-occurring variants
*across evidence types*. A somatic SNV co-occurring with an A-to-I edit on
one transcript cannot be represented by any haplotype panel — RNA editing
is not in the genome and never appears in a VCF of genotypes. ProHap
cannot see that proteoform in principle, not by omission.

The coverage rows are head-to-head measurements against these tools as run
here; the rest describes each tool's documented scope.

## What it does

- **Four evidence types in one database** — somatic SNVs/indels, A-to-I editing, gene fusions, alternative splicing.
- **Per-variant disposition** — every input gets a recorded outcome, including the ones that produce nothing.
- **[Combinatorial proteoforms](docs/USAGE.md#combinatorial-proteoforms)** — the protein carrying *all* co-occurring variants, phase honoured where the caller reports it.
- **[Non-canonical ORFs](docs/USAGE.md#about---include-noncanonical)** — three-frame translation of lncRNAs, pseudogenes and UTRs; opt-in, because it multiplies database size.
- **Nine release invariants** — mechanical self-checks before success is reported; any error-severity violation exits non-zero.
- **Provenance graph** — one record per sequence naming every variant behind it, and one row per novel peptide.
- **Byte-identical reproducibility** — two runs from one config produce the same release, verified over the whole tree.
- **Multi-assembly and multi-species** — GRCh38, GRCh37 and T2T-CHM13; human, mouse, rat and zebrafish.

## Measured

| check | result |
|---|---|
| ClinVar benchmark, 263 asserted rows | **99.2%** recall and precision |
| Ensembl VEP consequence agreement | **334/334** on shared transcripts |
| vs pypgatk 0.0.24, locus coverage | **260/260** against 259/260 |
| vs ProVar, variants represented | **2,865** against 911 |
| vs ProVar, entries carrying more than one variant | **7** against 0 |
| tests | **346**, offline, seconds |

```bash
make test && make reproducibility
```

Every disagreement, every caveat, and how to rerun each comparison:
[benchmarks/](benchmarks/).

## Limitations

- **Combining is not unique to v2p.** ProHap does it for phased germline
  haplotypes and does it well. What is unique here is combining *across*
  evidence types — see [why another one](#why-another-one).
- **Cross-evidence combination is not yet demonstrated on real data.** In
  HCC1395, 40 transcripts carry both a somatic variant and an editing
  site, but in none of them do both recode, so the measured count is 0.
  The capability is proven by test and by a
  [constructed case](benchmarks/cross_evidence_case.md), not by an
  observation.
- **Unphased combinations are hypotheses.** Where a caller reports no
  phase, co-occurrence is assumed, not known. Entries say which they are
  in `PHASE=`, and `--no-allow-unphased` drops them.
- **2 benchmark disagreements remain**, both the same FGFR3 variant
  (duplicated in the truth set) where isoform numbering differs.
- **Only genetic-code tables 1 and 2 are implemented.** A species needing
  the invertebrate (5) or yeast (3) mitochondrial table is **refused**
  rather than silently translated with the wrong one.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — installation, options, output, troubleshooting
- [docs/FORMAT_SPEC.md](docs/FORMAT_SPEC.md) — header grammar, all four styles
- [benchmarks/](benchmarks/) — truth sets, tool comparisons, how to rerun them
- [LICENSE](LICENSE) — MIT

# v2p — variant calls to protein sequences

Turns genomic and transcriptomic variant calls into a protein FASTA for
mass-spectrometry search, and records what happened to **every** input
variant — including the ones that produce nothing.

Reads four evidence types into one database: somatic SNVs/indels, A-to-I
RNA editing, gene fusions, and alternative splicing. No other tool covers
that combination.

## Install

```bash
pip install .
bash scripts/00_fetch_references.sh ref/     # GRCh38 + GENCODE v44, ~18 GB, once
```

Python 3.10–3.13. Only runtime dependency: `pyfaidx`.

## Use

```bash
v2p detect examples/                                    # what's in a folder
v2p run examples/ --ref ref/ --outdir out/ --name MYSAMPLE --logdir logs/
```

Files are identified by **content**, so no naming convention is needed.
`examples/` ships with the repo and works immediately.

Full guide: **[docs/USAGE.md](docs/USAGE.md)**.

## Output

| file | what it is |
|---|---|
| `<name>.target.fasta` | the database — search this |
| `<name>.target_decoy.fasta` | same, with decoys |
| `<name>.entries.tsv` | one row per sequence |
| `tables/provenance.jsonl` | one record per sequence, listing every variant that produced it |
| `tables/peptide_provenance.tsv` | one row per novel peptide, and whether it spans the variant residue |
| `_work/tables/disposition.tsv` | one row per **input** variant and its fate |
| `by_class/` | one FASTA per variant type |

## Where it sits

| tool | small variants | splice junctions | RNA editing | fusions | non-canonical ORFs |
|---|---|---|---|---|---|
| customProDB | yes | yes | no | no | no |
| QUILTS | yes | yes | no | no | partial |
| ProteoDisco | yes | yes | no | no | no |
| pypgatk | yes | no | no | no | yes |
| **v2p** | yes | yes | **yes** | **yes** | yes (opt-in) |

## Measured

| check | result |
|---|---|
| ClinVar benchmark, 263 asserted rows | **99.2% recall, 99.2% precision** |
| same, all 336 rows including unasserted | 99.4% recall, 99.4% precision |
| Ensembl VEP consequence agreement | **334/334 (100%)** on shared transcripts |
| A-to-I editing positive controls | **5/5** |
| GENCODE translation agreement | **100%** (single-transcript build) |
| vs pypgatk 0.0.24, locus coverage | **260/260** vs 259/260 |
| tests | **278**, offline, seconds |

```bash
make test                                   # 278 assertions, no reference needed
make reproducibility                        # two runs from one config, byte-identical
python3 scripts/09_benchmark.py --ref ref/  # the benchmark
v2p audit --ref ref/                        # check the reference interface
```

Nine release invariants run before `v2p run` reports success; any
error-severity violation exits non-zero.

## How it works

Variants apply to the **mature transcript**, not to genomic sequence —
that is what lets one code path serve DNA variants and RNA editing.
Translation runs from the annotated start to the first stop; frameshifts
read into the 3′ UTR.

Genomic coordinates are 1-based inclusive; transcript and CDS offsets are
0-based from the 5′ end. The genetic code is hard-coded in `seqops.py` so
a library update cannot change a translation silently. Transcript ties
break deterministically: MANE_Select → Ensembl_canonical → basic →
longest CDS → longest transcript → id.

## Limitations

- **Human and mouse only.** Other species need a `config/species/*.yaml`.
- **No ProteoDisco comparison.** It needs R, Bioconductor, BSgenome and a
  TxDb; not run. `benchmarks/compare_tools.md` has the procedure.
- **2 benchmark disagreements remain**, both the same FGFR3 variant
  (duplicated in the truth set) where the isoform numbering differs.
  Listed in `benchmarks/results/benchmark.md`.
- **73 of 336 ClinVar truth rows carry no assertion criteria.** They are
  excluded from the headline figure as conservatism about the truth data,
  not to flatter the result: scoring all 336 gives 99.4%/99.4%, slightly
  *higher* than the 263-row figure, with the same two disagreements.
  `--min-review-status any` reproduces it.
- **Non-canonical ORFs are off by default.** Enabling them takes the
  example database from 20,531 to 123,404 sequences; an inflated search
  space costs sensitivity at fixed FDR.

## Licence

MIT — see [LICENSE](LICENSE).

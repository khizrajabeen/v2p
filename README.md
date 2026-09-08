# v2p — variant calls to protein sequences

**v2p** builds mass-spectrometry search databases from variant calls, and
records what happened to **every** input variant — including the ones that
produce no protein. It reads DNA variants, RNA editing, gene fusions and
alternative splicing into one FASTA with a shared header vocabulary.

Silence is where wrong answers hide, so nothing is dropped without a
recorded reason.

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
- **Non-canonical ORFs** — three-frame translation of lncRNAs,
  pseudogenes and UTRs. Opt-in; see the warning below.
- **Multi-species** — human, mouse, rat and zebrafish ship; others are one
  YAML file.
- **Search-engine compatible** — four header conventions (UniProt, PEFF,
  pVAC, descriptive) and three decoy strategies.

## Installation

```bash
pip install .
```

Python 3.10–3.13. The only runtime dependency is `pyfaidx`.

Then fetch a reference once:

```bash
bash scripts/00_fetch_references.sh ref/      # GRCh38 + GENCODE v44, ~18 GB
```

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
| `--species <name\|yaml>` | human, mouse, rat, zebrafish, or a YAML path |

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

```
out/
  <name>.target.fasta              the database — search this
  <name>.target_decoy.fasta        the same, with decoys
  <name>.entries.tsv               one row per sequence
  <name>.summary_by_class.tsv      composition per variant type
  by_class/                        one FASTA per variant type
  tables/provenance.jsonl          every variant that produced each sequence
  tables/peptide_provenance.tsv    every novel peptide and what explains it
  qc/                              validation and recovery reports
  METHODS.md   MANIFEST.txt
  _work/tables/disposition.tsv     one row per INPUT variant and its fate
```

## Measured

| check | result |
|---|---|
| ClinVar benchmark, 263 asserted rows | **99.2% recall, 99.2% precision** |
| same, all 336 rows including unasserted | 99.4% / 99.4% |
| Ensembl VEP consequence agreement | **334/334 (100%)** on shared transcripts |
| A-to-I editing positive controls | **5/5** |
| GENCODE translation agreement | **100%** (single-transcript build) |
| vs pypgatk 0.0.24, locus coverage | **260/260** against 259/260 |
| tests | **287**, offline, seconds |

```bash
make test                                   # 287 assertions, no reference needed
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

## Limitations

- **Vertebrates only.** Human, mouse, rat and zebrafish ship; any other
  vertebrate is one `config/species/*.yaml`. Only genetic-code tables 1
  and 2 are implemented, so a species needing the invertebrate (5) or
  yeast (3) mitochondrial table is **refused** rather than silently
  translated with the wrong one.
- **Non-canonical ORFs are off by default**, and should stay off unless
  you are specifically hunting lncRNA or uORF peptides. Enabling them
  takes the example database from 20,531 to 175,552 sequences — 8.5×,
  after entries identical to a reference protein are already dropped. An
  inflated search space costs sensitivity at fixed FDR. A documented
  tradeoff, not a defect.
- **2 benchmark disagreements remain**, both the same FGFR3 variant
  (duplicated in the truth set) where isoform numbering differs.
- **No ProteoDisco number yet.** The driver and install route are in
  `benchmarks/compare_tools.md`; the run has not completed here.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — running your own data
- [docs/FORMAT_SPEC.md](docs/FORMAT_SPEC.md) — header grammar, all four styles
- [docs/TOOL_DESIGN.md](docs/TOOL_DESIGN.md) — input/output contract
- [benchmarks/](benchmarks/) — truth sets, results, tool comparison

## Licence

MIT — see [LICENSE](LICENSE).

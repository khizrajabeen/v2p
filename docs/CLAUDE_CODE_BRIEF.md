# Briefing for Claude Code

Paste the "Opening prompt" section into a fresh Claude Code session
started in the project root. The rest is reference material for that
session — point it at this file.

---

## Opening prompt

> You are taking over `v2p`, a working, validated bioinformatics pipeline
> that converts genomic and transcriptomic variant calls into amino-acid
> sequences. It is not a prototype: it has produced a delivered dataset,
> has 91 passing tests, and agrees with GENCODE's own translations at
> 100% on its single-transcript build.
>
> Read these three files before touching anything:
>   - `docs/BUILD_SPEC.md` — the seven-milestone plan and the competitive
>     landscape
>   - `docs/TOOL_DESIGN.md` — the input/output contract and architecture
>   - `README.md` — methods and known limitations
>
> Then run `python3 tests/test_pipeline.py` and confirm 91 passing.
>
> Your task is **Milestone 1 only**: implement `src/v2p/invariants.py` and
> wire it into `scripts/v2p run`, per the contract in BUILD_SPEC section
> "Milestone 1". Do not start Milestone 2 until M1's acceptance test
> passes.
>
> Hard constraints:
>   - Do not change translation behaviour. M1 is purely additive. Prove
>     it: run the pipeline before and after and diff the output FASTA.
>   - Every check you add gets both a positive and a negative test.
>   - No network access in tests. Fixtures are synthetic or committed.
>   - When you find a bug in existing code, fix it, add a regression test
>     named for the bug, and tell me — do not silently work around it.
>
> Ask before restructuring anything. The module boundaries were chosen
> deliberately and are documented.

---

## Where everything is

```
/mnt/d/conversion/
  HCC1395_variant_to_protein/          <- project root, run everything here
    src/v2p/                           the library
      discover.py      content-based input detection      [v2.2, new]
      annotation.py    GTF models, coordinate mapping     [core]
      seqops.py        genetic codes, translation         [core]
      nmd.py           nonsense-mediated decay            [core]
      peptides.py      digestion, novelty, decoys         [core]
      validate.py      UniProt/GENCODE comparison         [core]
      fasta.py         header styles, writers             [core]
      provenance.py    logging, checksums, environment    [core]
      parse/inputs.py  one adapter per evidence type
      build/           smallvar.py, fusion.py, splicing.py
    scripts/
      v2p                        the CLI: detect / run / audit
      00_audit_references.py     reference-interface audit
      00_fetch_references.sh     reference download + build assertion
      01_parse_inputs.py         -> unified variant manifest
      02_build_protein_fasta.py  translation, the core stage
      03_qc_report.py            input QC
      04_validate_uniprot.py     validation gates
      05_compare_tools.py        FASTA-vs-FASTA comparison
      06_package_release.py      packaging, decoys, manifest
      08_recovery_report.py      per-variant accounting
      clean.sh                   remove intermediates
    tests/
      make_fixture.py            synthetic genome generator
      test_pipeline.py           91 assertions
      fixtures/                  generated, not committed
    docs/
      BUILD_SPEC.md              THE PLAN — start here
      TOOL_DESIGN.md             input/output contract
      FORMAT_SPEC.md             header grammar, all four styles
      PLAN.md                    original conversion plan
    data/                        symlinks to the input files
    ref/                         GRCh38 + GENCODE v44 + UniProt (~18 GB)
    results/                     working state
    benchmark/                   a built database
    HCC1395_variant_proteins/    the delivered output
    logs/                        provenance JSONs
    benchmark.sh, finalize.sh    one-command build and rename
```

### Reference files (already downloaded, do not re-download)

```
ref/GRCh38.primary_assembly.genome.fa       3.1 GB, GRCh38 confirmed
ref/gencode.v44.annotation.gtf.gz           48 MB
ref/gencode.v44.pc_translations.fa.gz       12 MB, the validation ground truth
ref/uniprot_human_SP.fasta                  20,432 SwissProt entries
```

### Input files

```
data/sSNV_sIndel.vcf.gz                     33,221 somatic SNV/InDel, sites-only
data/HCC1395_high_confidence_RES_*.txt      8,093 A-to-I editing sites (ANNOVAR)
data/HCC1395_high_confidence_Fusion_*.csv   28 fusions
data/HCC1395_high_confidence_AS-LR_v1.csv   1,027 SUPPA2 splicing events
data/PGx_High-Confidence_Regions.bed        3.4M regions (unused so far)
```

---

## State of play

Working and validated:

- 42,369 input variants → 1,329 alter a protein → 2,150 distinct sequences
- GENCODE agreement **100%** (single-transcript build), 97.3% on the
  all-transcript build — the difference is `cds_end_NF` transcripts, where
  GENCODE truncates at the annotated CDS end and we read to the first stop
- nonsynonymous:synonymous ratio **2.75**, within the expected range
- ADAR positive controls recovered: GRIA2 Q607R, NEIL1 K242R, BLCAP Y2C,
  CDK13 Q103R, COG3 I635V
- reference audit reports **0 bugs**
- 91 tests, all offline, run in seconds

### Bugs found and fixed — the reason M1 matters

Every one of these produced *plausible-looking wrong output*, not an
error. None was caught by output-agreement testing.

| bug | effect | how it was found |
|---|---|---|
| transcripts looked up by gene symbol only | 33,221 sites-only VCF variants produced zero proteins | noticing "0 distinct genes" in a QC table |
| minus-strand AF/AL events dropped | 234 of 1,027 splicing events silently skipped | reading warning lines in a log |
| GTF frame column parsed, never used | 8,260 `cds_start_NF` transcripts translated out of frame | auditing the reference interface |
| selenocysteine annotation ignored | all 25 selenoproteins truncated at their first UGA | same audit |
| mitochondrial code not applied | chrM would read through its real stops | same audit |
| UTR variants classified as synonymous | inverted the nonsyn:syn ratio; 16 indels shipped as false "frameshifts" | the ratio looked biologically impossible |
| ambiguous gene symbols resolved by name | silent false negatives on paralogous loci | audit reported 25 ambiguous names |
| double-counted wild-types across files | per-type files did not sum to the combined file | the recipient asked why the numbers disagreed |

Read that table before writing any check. The pattern is consistent: the
output looked fine, and only an independent invariant or an audit of an
assumption exposed the error.

---

## Advanced features to build, in order

Full contracts are in `docs/BUILD_SPEC.md`. Summary of intent:

**M1 — output invariants.** Nine mechanical checks, wired into `run`,
exit non-zero on violation. Six of the eight bugs above are in this class.

**M2 — config and reproducibility.** YAML config sufficient to reproduce
a run. Acceptance test: two runs from one config produce byte-identical
FASTA. Catches set-iteration order, unsorted globs, embedded timestamps.

**M3 — ground-truth benchmark.** ProteoDisco published a 28-variant
COSMIC comparison; without an equivalent, no accuracy claim is
falsifiable. Include an ADAR recoding truth set — **no competing tool can
attempt it**, and that is the strongest single argument for this tool
existing.

**M4 — species independence.** Remove hard-coded `_HUMAN`, `OX=9606`,
GRCh38 contig lengths. A pure refactor: prove human output is
byte-identical afterwards.

**M5 — non-canonical ORFs.** Three-frame translation of lncRNAs and
pseudogenes. The one capability where pypgatk is genuinely ahead. Off by
default: it multiplies database size, and an inflated search space costs
sensitivity.

**M6 — distribution.** `pip install v2p`, bioconda recipe, GitHub Actions
across Python 3.10–3.13, Dockerfile pinning the reference release.

**M7 — the differentiator.** Provenance graph (which calls produced this
sequence), peptide-level traceability (which variants explain this
peptide), continuous reference auditing in CI, and an optional VEP
consequence cross-check. This is the part worth publishing.

---

## Things to know before changing code

- **Coordinates.** Genomic are 1-based inclusive; transcript and CDS are
  0-based offsets from the 5' end in transcription direction. Every
  conversion goes through `annotation.py` so the convention lives in one
  place. BED is 0-based half-open — never mix it in by hand.
- **Variants apply to the mature transcript, not to genomic sequence.**
  That is what lets one code path serve DNA variants and RNA editing.
  Do not "simplify" this into genomic-space substitution.
- **The genetic code is hard-coded in `seqops.py`, deliberately.** A
  library update must not be able to change a translation silently.
- **Determinism is a requirement.** Transcript selection ties break on
  MANE_Select → Ensembl_canonical → basic → longest CDS → longest
  transcript → lexicographic id. Never leave a tie unbroken.
- **GENCODE agreement is the correctness gate, not UniProt.** UniProt and
  GENCODE choose canonical isoforms independently and often disagree on
  the start codon; that is an annotation difference, not an error.

## Do not do these

- Do not add MHC binding prediction or neoantigen ranking — pVACtools
  owns that and does it well.
- Do not try to match VEP's consequence catalogue.
- Do not enable non-canonical ORFs by default.
- Do not remove the audit. If a reference-interface assumption changes,
  the audit must fail before the output does.
- Do not report only favourable benchmark results. If ProteoDisco beats
  us on small variants, that belongs in the README.

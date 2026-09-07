# v2p — build specification

Handoff document. Written to be executed by an agent with repo access.
Each milestone has a contract, an acceptance test, and a definition of
done. Do not proceed to the next milestone until the current one's
acceptance test passes.

---

## 0. Competitive position — read this first

Four tools occupy this space. Their capabilities, verified from the
publications:

| tool | year | language | small variants | splice junctions | RNA editing | fusions | non-canonical ORFs | workflow | status |
|---|---|---|---|---|---|---|---|---|---|
| customProDB | 2013 | R | yes | yes | no | no | no | none | unmaintained |
| QUILTS | 2016 | Python | yes | yes | no | no | partial | none | dormant |
| ProteoDisco | 2021 | R/Bioconductor | yes (SNV/MNV/InDel) | yes | no | no | no | none | maintained |
| pypgatk / pgdb | 2021 | Python | yes | no | no | no | **yes** (3-frame) | Nextflow, bioconda | maintained |
| **v2p** | — | Python | yes | yes | **yes** | **yes** | not yet | none yet | working |

### What we already do that none of them do

1. **RNA editing as a first-class evidence type.** No tool above accepts
   an editing table. This matters: A-to-I recoding produces proteoforms
   with no genomic basis, so they are invisible to any DNA-driven
   pipeline. Our GRIA2 Q607R and NEIL1 K242R recoveries are the proof.
2. **Gene fusions.** AGFusion does fusions, but it does not produce a
   search database, and no database tool ingests fusion calls.
3. **Per-variant disposition.** Every input variant gets a recorded
   outcome. The others emit what they produced and are silent about what
   they did not. Silence is where wrong answers hide.
4. **Reference-interface auditing.** Verifying the GTF, genome and
   proteome are being *read* correctly, separately from checking that the
   output looks right. This found four real bugs — ignored GTF frame,
   selenocysteine truncation, mitochondrial code, ambiguous gene symbols
   — that output-agreement testing missed entirely.
5. **Content-based input detection.** Point it at a folder.

### What they have that we lack, and must acquire

| gap | who has it | milestone |
|---|---|---|
| non-canonical / cryptic ORFs (3-frame translation of lncRNA, pseudogenes) | pypgatk | M5 |
| published ground-truth benchmark (ProteoDisco used 28 COSMIC variants) | ProteoDisco | M3 |
| installable via pip and conda | pypgatk | M6 |
| workflow-manager integration | pgdb (Nextflow) | M6, optional |
| multi-species | pypgatk | M4 |
| peer-reviewed publication | all | after M6 |

**The claim to make, and only this one:** the first tool to build a single
protein database from DNA variants, RNA editing, fusions and splicing
together, with every input variant's fate recorded. Not "better than
VEP" — VEP is not a database builder. Not "most accurate" — that needs
the M3 benchmark to say.

---

## Milestone 1 — output invariants

**Why first.** Six real errors reached a shipped file during development:
double-counted wild-types, per-type files not summing to the combined
file, UTR variants labelled frameshift, an incomplete reference proteome,
duplicated accessions, and a stale manifest. Every one is mechanically
checkable. A tool that can ship these is not trustworthy regardless of
what else it does.

**Contract.** New module `src/v2p/invariants.py`:

```python
@dataclass
class Violation:
    check: str
    severity: str        # "error" | "warning"
    message: str
    examples: list[str]  # at most 5 offending ids

def check_release(outdir: Path, manifest: Path | None,
                  disposition: Path | None) -> list[Violation]: ...
```

Checks, all mandatory:

| id | check | severity |
|---|---|---|
| I1 | per-type files partition the combined file: union equals it, pairwise intersections empty | error |
| I2 | every input variant appears exactly once in the disposition | error |
| I3 | variant + reference + wildtype sequence counts sum to the combined file | error |
| I4 | every sequence id is unique | error |
| I5 | sequence line layout matches the supplied proteome template | error |
| I6 | reference translations agree with the annotation source above threshold | error |
| I7 | no entry labelled as a sequence-changing consequence has a sequence equal to its own reference | warning |
| I8 | `MANIFEST.txt` checksums match the files on disk | error |
| I9 | decoy count equals target count when decoys were requested | error |

**Wiring.** `v2p run` calls `check_release` before printing success. Any
error-severity violation exits non-zero and prints the offending ids.

**Acceptance test** `tests/test_invariants.py`: construct a miniature
release, assert zero violations; then inject each of the nine faults in
turn and assert exactly the corresponding check fires. Nine positive and
nine negative cases.

**Done when** the suite passes and `v2p run` on the HCC1395 data reports
zero violations.

---

## Milestone 2 — configuration and reproducibility

**Contract.** `src/v2p/config.py` loading the YAML schema in
`docs/TOOL_DESIGN.md`, with defaults, type validation and clear errors on
unknown keys. `v2p run --config run.yaml` must be sufficient to reproduce
a run; `v2p run --write-config` emits the config a command line implies.

**Reproducibility test.** Run twice into different directories from the
same config; assert the two combined FASTA files are byte-identical.
This is the test that catches accidental nondeterminism — set iteration
order, unsorted globs, timestamps in sequence files.

**Done when** the double-run test passes and every config key appears in
the provenance JSON.

---

## Milestone 3 — the ground-truth benchmark

**Why this matters most for credibility.** ProteoDisco's paper compares
against customProDB and QUILTS on 28 COSMIC variants with known protein
consequences. Without an equivalent, "more accurate" is unfalsifiable.

**Contract.** `benchmarks/` holding:

1. `truth/cosmic_variants.tsv` — 50 or more variants with an
   independently established protein consequence, drawn from COSMIC or
   ClinVar. Columns: chrom, pos, ref, alt, gene, transcript, expected
   HGVS protein change, source, source version.
2. `truth/editing_sites.tsv` — the canonical ADAR recoding sites
   (GRIA2 Q607R, NEIL1 K242R, BLCAP Y2C, CDK13 Q103R, COG3 I635V, and
   others from REDIportal), with expected consequences. **No other tool
   can attempt this set.**
3. `truth/fusions.tsv` — characterised fusions with published breakpoints
   and protein consequences (BCR-ABL1, EML4-ALK, TMPRSS2-ERG).
4. `scripts/09_benchmark.py` — runs the pipeline over the truth set and
   reports per-category precision and recall, plus every disagreement
   with its expected and observed value.
5. `benchmarks/compare_tools.md` — same truth set through ProteoDisco and
   pypgatk where they support the category, with install commands and
   exact versions so the comparison is reproducible.

**Report honestly.** If a competitor matches or beats us on small
variants, publish that. The differentiator is coverage of editing and
fusions, not superiority on the shared subset. A benchmark that only
flatters is worthless.

**Done when** `python scripts/09_benchmark.py` prints a per-category
table and the small-variant recall is at or above 95%.

---

## Milestone 4 — species and annotation independence

Currently hard-coded: `_HUMAN` in entry names, `OS=Homo sapiens OX=9606`,
GRCh38 contig lengths in the audit, GENCODE-specific tag names.

**Contract.** A `Species` object carrying scientific name, taxon id,
entry-name suffix, expected contig lengths and mitochondrial code table,
loaded from `config/species/*.yaml`. Ship human and mouse; the audit must
degrade gracefully for a species with no reference lengths.

**Acceptance test.** Build a mouse fixture (synthetic, same construction
as the human one) and assert correct headers, `OX=10090`, and that the
build assertion does not falsely fail.

**Done when** the human output is byte-identical to before the change —
this is a refactor and must prove it changed nothing.

---

## Milestone 5 — non-canonical ORFs

The one capability gap where pypgatk is genuinely ahead: three-frame
translation of lncRNAs, pseudogenes and untranslated regions, to catch
proteins from "dormant" genomic regions.

**Contract.** `src/v2p/build/noncanonical.py`:

- three-frame translation of transcripts whose biotype is non-coding
- configurable minimum ORF length, default 30 residues
- ORF start policy: first ATG, or any-start with the frame recorded
- new variant classes `NC_LNCRNA`, `NC_PSEUDOGENE`, `NC_UTR`, `NC_ALTORF`
- these entries must be **excluded by default**. They can multiply
  database size several-fold, and an inflated search space costs
  sensitivity. `--include-noncanonical` opts in.

**Acceptance test.** A fixture lncRNA with a known 40-residue ORF in
frame 2; assert it is found, that frame is recorded, that a 20-residue
ORF is rejected at the default threshold, and that the entries are absent
unless the flag is given.

**Done when** the flag works and the default output is unchanged.

---

## Milestone 6 — distribution

**Contract.**

- `pyproject.toml`, console entry point `v2p = v2p.cli:main`
- `pip install v2p` works in a clean virtualenv
- conda recipe targeting bioconda
- GitHub Actions: test suite on Python 3.10 through 3.13, no reference
  download, under two minutes
- `Dockerfile` pinning the reference release
- optional `nextflow/main.nf` wrapping `v2p run`, for parity with pgdb

**Done when** `pip install .` in a clean environment followed by
`v2p detect <folder>` works with no `PYTHONPATH` manipulation.

---

## Milestone 7 — the differentiator, built out

Everything above is parity work. This is the part no other tool has.

### 7a. Provenance graph

Currently the variant-to-sequence mapping is recoverable by joining two
TSVs. Make it first class: `tables/provenance.jsonl`, one record per
output sequence, listing every input variant contributing to it, the
transcript, and the reference sequence it derives from. A user should be
able to ask "which calls produced this peptide" and get an answer without
writing a join.

### 7b. Peptide-level provenance

For each novel tryptic peptide: which sequences contain it, which
variants explain it, and whether it spans a variant residue. This is what
a proteomics group actually needs after a search returns a hit, and no
current tool provides it.

### 7c. Continuous reference auditing

Run `v2p audit` in CI against the pinned reference. When GENCODE
releases a new version, the audit tells you what changed about the
*interface* — new feature types, changed tag vocabulary, altered frame
conventions — before it silently changes your output.

### 7d. Cross-tool consequence check

When VEP is installed, compare consequence assignments on the shared
subset and report disagreements. Not as validation — as a second opinion
worth naming in a methods section.

---

## Order of work

M1 and M2 first: they make the tool trustworthy and reproducible, and
both are small. M3 next, because without a benchmark no accuracy claim is
defensible. M4 and M6 are mechanical. M5 only if lncRNA-derived peptides
are actually wanted — it is the largest single addition and it degrades
default output if enabled carelessly. M7 last, and it is the part worth
publishing.

## Constraints for the implementing agent

- **Do not change translation behaviour** while doing M1, M2, M4 or M6.
  Those are refactors. Prove it with a byte-identical output check.
- **Every bug fixed gets a regression test**, named for the bug.
- **No network access in tests.** All fixtures synthetic or committed.
- **The audit is not optional.** If a reference-interface assumption
  changes, the audit must fail before the output does.
- **Report negative results.** If a benchmark shows a competitor ahead,
  that goes in the README, not a drawer.

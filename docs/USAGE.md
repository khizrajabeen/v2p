# Running your own dataset

Everything the README does not have room for: installation, the full
option reference, what each output file holds, and what to do when a run
stops.

- [Installation](#installation)
- [The short version](#the-short-version)
- [Step 1 — put the reference in place, once](#step-1--put-the-reference-in-place-once)
- [Step 2 — put your calls in a folder](#step-2--put-your-calls-in-a-folder)
- [Step 3 — check what was detected](#step-3--check-what-was-detected)
- [Step 4 — convert](#step-4--convert)
- [What you get](#what-you-get)
- [Reading a header](#reading-a-header)
- [Command reference](#command-reference)
  - [Main interface](#main-interface)
  - [Choosing a reference](#choosing-a-reference)
  - [Naming inputs explicitly](#naming-inputs-explicitly)
  - [Output control](#output-control)
  - [Reproducibility](#reproducibility)
  - [Flags worth knowing](#flags-worth-knowing)
- [Combinatorial proteoforms](#combinatorial-proteoforms)
- [About `--include-noncanonical`](#about---include-noncanonical)
- [Reproducing a run](#reproducing-a-run)
- [Non-human data](#non-human-data)
- [How it works](#how-it-works)
- [Analysis scripts](#analysis-scripts)
- [When something goes wrong](#when-something-goes-wrong)
- [What the tool refuses to do](#what-the-tool-refuses-to-do)
- [Deliberately out of scope](#deliberately-out-of-scope)
- [Checking the reference itself](#checking-the-reference-itself)

## Installation

```bash
pip install .
```

Python 3.10–3.13. The only runtime dependency is `pyfaidx`.

## The short version

```bash
v2p detect  /path/to/your/calls               # what is in there
v2p run     /path/to/your/calls --ref ref/ \
            --outdir out/ --name MYSAMPLE --logdir logs/
```

`detect` writes nothing, so it is safe to point at anyone's data to see
what v2p would make of it. Always run it first — it says which files were
recognised, how confident it is, and what it would do.

## Step 1 — put the reference in place, once

```bash
bash scripts/00_fetch_references.sh ref/                     # GRCh38, the default
bash scripts/00_fetch_references.sh ref37/  --assembly GRCh37
bash scripts/00_fetch_references.sh reft2t/ --assembly T2T
```

Each downloads the genome, annotation and — where the source publishes
them — the annotation's own protein translations, then asserts the build
by contig length before you spend an hour translating against the wrong
one. `GENCODE_RELEASE=45 bash scripts/00_fetch_references.sh ref/` picks a
different release. Do it once; it is the same for every sample.

You also need a reference proteome FASTA (UniProt SwissProt for your
species) in the same directory. `detect --ref ref/` will say if it is
missing.

## Step 2 — put your calls in a folder

**No naming convention is required.** Files are identified by content, so
call them whatever your collaborator called them. Drop any mix of these
into one directory:

| what you have | recognised by | required fields |
|---|---|---|
| somatic SNVs / indels | a `##fileformat=VCF` header | CHROM, POS, REF, ALT |
| A-to-I RNA editing | tab-separated with `Func.refGene` or `Func.ensGene` | Chr, Start, Ref, Alt |
| gene fusions | delimited text with two gene columns and two `chr:pos` columns | gene1, gene2, breakpoint1, breakpoint2 |
| alternative splicing | a column of SUPPA2 event ids (`<gene>;<TYPE>:<chr>:...:<strand>`) | event id, optionally dPSI |
| confidence regions | BED, three or more columns | chrom, start, end |

You do not need all four. One VCF on its own is a perfectly good input.
See `examples/` for a working folder.

## Step 3 — check what was detected

```bash
v2p detect mycalls/ --ref ref/
```

```
file                       identified as                conf  records
somatic.vcf                somatic SNV / InDel calls     100%   33,221
editing.txt                RNA editing sites              98%    8,093
fusions.csv                gene fusion calls              95%        28

run plan
  somatic SNV / InDel calls              somatic.vcf
  reference genome                       GRCh38.primary_assembly.genome.fa
  genome build                           GRCh38

no problems. `v2p run` would convert somatic SNV / InDel calls, ...
```

If two files are equally plausible for one role it says so and refuses to
guess — name the one you want rather than letting it pick. `-v` shows the
evidence behind every decision.

## Step 4 — convert

```bash
v2p run mycalls/ --ref ref/ --outdir out/ --name MYSAMPLE --logdir logs/
```

A couple of minutes for a whole-exome call set. It prints the five stages,
then the invariant report, and **exits non-zero if the release contradicts
itself**.

Keep `--logdir` outside `--outdir`. `MANIFEST.txt` is written last, so
anything added to the release afterwards makes it stale and the I8
invariant will correctly fail the run.

## What you get

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

`disposition.tsv` is the one people forget. Every input variant appears in
it exactly once, including the ones that produced nothing, with the reason
— intronic, no coding transcript, synonymous, and so on. It is the only
way to state a recovery rate honestly.

## Reading a header

```
>vr|MYSAMPLE_SNV_000001|AGRN_HUMAN_SNV Agrin (p.E1608Q; somatic
 single-nucleotide variant) OS=Homo sapiens OX=9606 GN=AGRN VT=SNV
 CSQ=missense TX=ENST00000379370.7 PC=p.E1608Q LOC=chr1:1049980G>C
 POS=1608 NMD=escapes NOVELPEP=16 VARPEP=4 NPEP=403
```

`vr|` is a variant entry, `sp|` a reference one, so a search engine can
separate them on the database field alone. `VT=` variant type, `CSQ=`
consequence, `PC=` protein change, `NOVELPEP=` peptides absent from the
reference proteome, `VARPEP=` how many of those span the variant residue.
Combinatorial entries also carry `PHASE=phased` or `PHASE=unphased`.

Full grammar in [FORMAT_SPEC.md](FORMAT_SPEC.md).

## Command reference

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
| `--combine-max <n>` | most variants combined on one transcript, default 8 |
| `--no-allow-unphased` | keep only combinations the caller actually phased |
| `--min-af <f>` | drop variants below this INFO allele frequency |
| `--max-combinatorial-fraction <f>` | fail rather than inflate the database past this share |
| `--include-noncanonical` | three-frame translate non-coding transcripts |
| `--split-by-type` / `--no-split` | per-variant-type FASTA files |

### Reproducibility

| option | description |
|---|---|
| `--config <yaml>` | run from a config file |
| `--write-config <yaml>` | emit the config a command line implies |
| `--min-agreement <f>` | translation-agreement gate, default 0.90 |
| `--skip-invariant <ID>` | turn off one release check, recorded in provenance |

### Flags worth knowing

| flag | why |
|---|---|
| `--transcript-mode representative` | one transcript per gene instead of all. Much smaller database. |
| `--header-style peff \| uniprot \| pvac \| descriptive` | re-emit in another convention without re-translating |
| `--decoys none` | if your search engine generates its own |
| `--drop-unchanged` | drop synonymous and UTR variants (kept by default) |
| `--species mouse` | non-human. See `config/species/`. |
| `--include-noncanonical` | three-frame translate lncRNAs and pseudogenes. **Read the warning.** |
| `--min-agreement 0.95` | raise the translation-agreement gate |
| `--skip-invariant I5` | turn off one release check, recorded in the provenance JSON |

## Combinatorial proteoforms

`--combine-variants` emits, for each transcript carrying two or more
co-occurring variants, the protein carrying *all* of them.

The reason is peptide-level. A tryptic peptide spanning two variants is in
neither single-variant entry nor the reference, so a database built one
variant at a time cannot identify it at any FDR.
[ProHap](https://doi.org/10.1038/s41592-024-02506-0) measured this at
12.4% of substitutions for common germline haplotypes; on HCC1395's
sparser somatic calls it is 2.9%.

Two variants can even share a codon. In HCC1395, FKTN carries
`p.[D225K]` — a lysine neither single-variant entry produces, and one
that creates a tryptic cleavage site, changing the peptides on both
sides of it.

Three rules keep this from inflating the database:

- **Every entry must earn its place.** A combination is emitted only if
  at least one of its tryptic peptides is absent from the reference *and*
  from every single-variant protein of the same transcript. On HCC1395
  that drops 40 of 49 candidates and keeps 9.
- **Phase is honoured where the caller reports it.** With `GT` and `PS`,
  variants combine within a haplotype, and two on opposite haplotypes are
  never combined — that combination is not uncertain, it is false. A
  homozygous variant joins both haplotypes.
- **Without phase, the entry is a hypothesis**, marked `PHASE=unphased`
  in its header so a search can filter it out later.
  `--no-allow-unphased` drops them at build time instead.

Protein changes use HGVS allele notation, which distinguishes the two
cases: `p.[A12V;G45S]` when phased, `p.[A12V(;)G45S]` when not.

`--combine-max` caps how many variants are combined on one transcript
(default 8); above that, a transcript is far likelier to be an alignment
artefact than a real proteoform.

## About `--include-noncanonical`

Three-frame translation of lncRNAs, pseudogenes, and the 5′/3′ UTRs of
coding transcripts (uORFs and dORFs), with the reading frame recorded per
entry and a configurable length floor. Entries identical to a reference
protein are dropped rather than duplicating the database.

It is off by default for a reason. On the example data it takes the
database from 20,531 to 123,404 sequences — six times larger. An inflated
search space costs sensitivity: more candidate peptides means a higher
score threshold at the same FDR, so you lose real identifications
elsewhere. Turn it on when you are specifically hunting lncRNA-derived
peptides, not by default.

## Reproducing a run

```bash
v2p run mycalls/ --ref ref/ --write-config run.yaml   # capture
v2p run --config run.yaml --outdir out2/              # replay
```

Two runs from one config produce a byte-identical release; `make
reproducibility` is that test. Every parameter also lands in
`logs/*.provenance.json` with input checksums, so a config file plus those
checksums fully determine the output.

## Non-human data

```bash
v2p run mycalls/ --ref ref_mouse/ --species mouse --outdir out/
```

Human, mouse, rat and zebrafish ship, along with GRCh37 and T2T profiles
for human. For anything else, drop a YAML into `config/species/` — copy
`config/species/mouse.yaml` — and pass its name or path. Only headers and
the audit's contig-length assertions change; the translation logic is
species-neutral.

A species needing a genetic-code table other than 1 (standard) or 2
(vertebrate mitochondrial) is **refused** rather than translated with the
wrong one.

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

## Analysis scripts

| script | description |
|---|---|
| `scripts/09_benchmark.py` | score against the truth sets, per category |
| `scripts/10_vep_annotate.py` | annotate a VCF via Ensembl VEP REST (no local cache) |
| `benchmarks/run_provar.sh` | run ProVar (ProHap's companion) on the same VCF |
| `benchmarks/compare_provar.py` | score v2p against a ProVar database |
| `benchmarks/cooccurrence_stat.py` | how often two substitutions share a peptide |
| `benchmarks/cross_evidence_case.py` | the worked cross-evidence proteoform |
| `benchmarks/run_proteodisco.R` | run ProteoDisco on the same reference |

## When something goes wrong

**"cannot run - missing: reference genome"** — the genome, annotation or
proteome was not found. Pass `--ref` pointing at the directory holding
them.

**"refusing to run; re-run with --force"** — detection found a problem,
usually two equally plausible files for one role. Read it, then either name
the file explicitly or pass `--force` if it is harmless.

**"N invariant violation(s)"** — the release contradicts itself and the run
exits non-zero. The report names the failing check. This is the tool
refusing to hand you a database it cannot vouch for; `--force` deliberately
cannot override these. `--skip-invariant <ID>` can, per check, and records
the fact.

**"validation below threshold"** — reference translations disagreed with
the annotation source's own by more than `--min-agreement` allows. Usually
a mismatched genome and GTF release.

## What the tool refuses to do

Some situations end the run rather than produce a database that looks
fine and is not:

- **mixed genome builds** — a GTF and genome from different assemblies
- **an undeterminable annotation release**, when the run is asked to be
  reproducible: a release that cannot be named cannot be recorded
- **an evidence file that parses but yields zero records** — silence
  there almost always means a format mismatch, not an empty call set

## Deliberately out of scope

MHC binding prediction, expression filtering, clinical interpretation,
and a graphical interface. Each is a different product, and pVACtools
already does the first two well.

## Checking the reference itself

```bash
v2p audit --ref ref/
```

Verifies that the GTF, genome and proteome are being *read* correctly,
separately from whether the output looks right. Four real bugs were found
this way that output-agreement testing missed entirely. Run it whenever you
change reference release.

# v2p — design specification

## What this is for

Convert genomic and transcriptomic variant calls into amino-acid sequences,
with every output sequence traceable to the call that produced it.

## Honest positioning

Claiming to beat every existing tool would be a marketing statement, not a
design goal. Here is what the landscape actually looks like:

| tool | does better than us | we should not try to match |
|---|---|---|
| Ensembl VEP | annotation breadth, plugin ecosystem, 15 years of edge cases, regulatory and non-coding consequence | its consequence catalogue |
| AGFusion | fusion transcript modelling, domain-level annotation of the chimera | its fusion domain analysis |
| pVACtools | MHC binding prediction, neoantigen ranking, clinical workflow | anything immunogenicity-related |
| customProDB | nothing current — unmaintained, Ensembl-89-era | — |

What none of them does, and what this tool should be built around:

1. **Four evidence types in one database.** VEP handles small variants,
   AGFusion handles fusions, splicing needs a separate pipeline, and RNA
   editing needs manual VCF conversion. Producing one coherent database
   from all four, with a shared header vocabulary, is genuinely unserved.

2. **Per-variant disposition.** Every input variant gets a recorded
   outcome — protein built, synonymous, UTR, intronic, intergenic, no
   transcript, error. Most tools emit what they produced and are silent
   about what they did not. Silence is where wrong answers hide.

3. **Reference-interface auditing.** Checking that the GTF, genome and
   proteome are being *read* correctly, separately from checking that the
   output looks right. This found four real bugs in our own pipeline that
   output-agreement testing missed entirely.

4. **Peptide-level accounting.** Reporting how many tryptic peptides an
   entry contributes that are absent from the reference, and how many of
   those actually span the variant residue.

Those four are the product. Everything else should be adequate, not
ambitious.

---

## Input contract

The tool accepts a **folder** and works out what is in it. No filename
conventions required — files are identified by content, because the same
data arrives named differently from every collaborator.

### Recognised evidence types

| type | detected by | required fields |
|---|---|---|
| small variants | `##fileformat=VCF` header | CHROM, POS, REF, ALT |
| RNA editing | tab-separated with `Func.refGene` or `Func.ensGene` | Chr, Start, Ref, Alt |
| gene fusion | delimited text with two gene columns and two `chr:pos` columns | gene1, gene2, breakpoint1, breakpoint2 |
| alternative splicing | column matching `<gene>;<TYPE>:<chr>:...:<strand>` | event id, optionally dPSI |
| confidence regions | BED, three or more columns | chrom, start, end |
| reference genome | FASTA, nucleotide alphabet | — |
| annotation | GTF/GFF with `transcript_id` | — |
| reference proteome | FASTA, protein alphabet, `sp\|`/`tr\|` headers | — |

Detection reports a confidence and the evidence for each decision, and
refuses to guess when two types are equally plausible. A run plan is
written before anything is processed, so the user can correct it.

### Input validation, before any work

- genome build inferred from VCF contig lengths and confirmed against the
  genome FASTA; mismatch is fatal
- contig naming (`chr1` vs `1`) checked for consistency across all files
- coordinate base recorded per format (VCF 1-based inclusive, BED 0-based
  half-open) and converted once, centrally
- reference alleles verified against the genome on a sample; a high
  mismatch rate aborts rather than producing silently wrong codons

### What the tool must refuse

- mixed genome builds
- a GTF whose release cannot be determined, when the caller asks for
  reproducible output
- silently proceeding when an evidence file parses but yields zero records

---

## Output contract

```
<outdir>/
  sequences/
    all_variant_types.fasta          every sequence, one file
    variant_sequences_only.fasta     variant-derived only
    reference_proteome.fasta         supplied proteome, verbatim
    wildtype_counterparts.fasta      wild-type of affected transcripts
    by_variant_type/                 one file per type + INDEX.txt
  tables/
    entry_details.tsv                one row per output sequence
    disposition.tsv                  one row per input variant
    summary_by_variant_type.tsv
    recovery_by_class.tsv
  reports/
    RECOVERY.txt   VALIDATION.txt   INPUT_SUMMARY.txt   AUDIT.txt
  logs/
    <stage>.<run_id>.log
    <stage>.<run_id>.provenance.json
  MANIFEST.txt
  README.txt
```

**Invariants the tool checks before declaring success:**

1. the per-type files partition the combined file exactly — no overlap, no
   gap
2. every input variant appears exactly once in `disposition.tsv`
3. `variant + reference + wildtype` counts sum to the combined file
4. every sequence id is unique
5. sequence layout matches the supplied proteome template (wrapped or not)
6. reference translations agree with the annotation source's own
   translations above a configurable threshold

A run that violates any invariant exits non-zero. These are the checks
that caught the double-counted wild-types and the mislabelled UTR
variants; they belong in the tool, not in a reviewer's head.

### Header vocabulary

```
>vr|<id>|<GENE>_HUMAN_<TYPE> <description> OS=... OX=... GN=...
   VT=<variant type>        SNV, INDEL, RNA_EDITING, FUSION, AS_*
   CSQ=<consequence>        missense, frameshift, 5_prime_UTR, ...
   TX=<transcript>  PC=<protein change>  LOC=<locus>  POS=<residue>
   NMD=<likely|escapes>     nonsense-mediated decay prediction
   NOVELPEP= VARPEP= NPEP=  peptide-level novelty counts
   ALSO=<types>             other classes giving the same sequence
```

---

## Command-line interface

Three verbs, not nine scripts:

```
v2p detect  <folder>                  what is in here, and what would run
v2p run     <folder> --ref <dir>      the whole conversion
v2p audit   --ref <dir>               check the reference interface only
```

`run` accepts a config file so a run is reproducible from one artefact:

```yaml
inputs:  auto                 # or explicit paths per evidence type
reference:
  genome: ref/GRCh38.fa
  annotation: ref/gencode.v44.gtf.gz
  proteome: ref/uniprot_human_SP.fasta
translate:
  transcript_mode: all        # all | representative
  keep_unchanged: true        # synonymous and UTR variants
  genetic_code: auto          # per-contig; table 2 on chrM
output:
  header_style: uniprot       # uniprot | peff | pvac | descriptive
  decoys: pseudo_reverse      # none | reverse | pseudo_reverse | shuffle
  split_by_type: true
validate:
  min_translation_agreement: 0.90
```

Every parameter lands in the provenance JSON, so a config file plus input
checksums fully determine the output.

---

## Architecture

```
v2p/
  discover.py     content-based input detection -> run plan
  config.py       config loading, defaults, validation
  cli.py          detect / run / audit
  annotation.py   transcript models, coordinate mapping   [exists]
  seqops.py       genetic codes, translation              [exists]
  nmd.py          decay prediction                        [exists]
  peptides.py     digestion, novelty, decoys              [exists]
  validate.py     UniProt/GENCODE comparison              [exists]
  fasta.py        header styles, writers                  [exists]
  invariants.py   the nine output checks                  [exists]
  parse/          one adapter per evidence type           [exists]
  build/          one builder per variant class           [exists]
```

Most of it exists and is tested. The work is the front end — detection,
config, one CLI — and the invariant checks.

### Adding a new evidence type

One adapter in `parse/` yielding the common record schema, one builder in
`build/` if the biology differs from existing classes, one entry in the
type vocabulary, and fixtures. Nothing else changes. That extensibility is
the point: circular RNA, retained-intron neoepitopes and structural
variants should each be an adapter, not a fork.

---

## Test strategy

Four layers, all offline:

1. **Unit** — genetic codes, coordinate conversion, digestion, decoys, NMD.
2. **Known-answer** — a synthetic genome where the expected protein for
   every case is known by construction, because the transcript is built
   first and then placed into the contigs. Plus and minus strand,
   selenoprotein, mitochondrial, duplicated gene symbol, `cds_start_NF`.
3. **Invariant** — the nine output checks, on a miniature end-to-end run,
   each with a positive case and an injected fault.
4. **Regression** — one test per bug ever found. Currently: minus-strand
   AF/AL ordering, GTF frame, selenocysteine, UTR-vs-synonymous,
   ambiguous symbols, positional lookup.

Nothing requires a reference download, so the suite runs in seconds in CI.

---

## Roadmap

| phase | work | why |
|---|---|---|
| 1 | `discover.py`, `config.py`, `cli.py` | any input folder, one command |
| 2 | `invariants.py` wired into `run` | the checks that caught real errors |
| 3 | consequence cross-check against VEP when installed | independent second opinion, optional |
| 4 | circular RNA and structural-variant adapters | the obvious next evidence types |
| 5 | packaging: `pip install v2p`, conda recipe, CI | so others can use it |

Phase 1 and 2 are the difference between a pipeline that works and a tool
someone else can run. Phases 3 to 5 are what makes it citable.

## Deliberately out of scope

MHC binding prediction, expression filtering, clinical interpretation, and
a graphical interface. Each is a different product, and pVACtools already
does the first two well.

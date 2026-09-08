# Running your own dataset

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
bash scripts/00_fetch_references.sh ref/
```

Downloads GRCh38 plus GENCODE v44 (~18 GB) and asserts the build. Do it
once; it is the same for every sample.

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

```
out/
  MYSAMPLE.target.fasta          the database — search this
  MYSAMPLE.target_decoy.fasta    the same, with decoys appended
  MYSAMPLE.entries.tsv           one row per sequence: class, gene,
                                 transcript, protein change, peptide counts
  MYSAMPLE.summary_by_class.tsv  composition at a glance
  by_class/                      one FASTA per variant type
  qc/                            validation and recovery reports
  METHODS.md                     what was done, in prose
  MANIFEST.txt                   SHA-256 of every file
  _work/tables/disposition.tsv   one row per INPUT variant and its fate
```

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

Full grammar in `docs/FORMAT_SPEC.md`.

## Flags worth knowing

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

### About `--include-noncanonical`

On the example data it takes the database from 20,531 to 123,404
sequences — six times larger. An inflated search space costs sensitivity:
more candidate peptides means a higher score threshold at the same FDR, so
you lose real identifications elsewhere. Turn it on when you are
specifically hunting lncRNA-derived peptides, not by default.

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

Ships with human and mouse. For anything else, drop a YAML into
`config/species/` — copy `config/species/mouse.yaml` — and pass its name or
path. Only headers and the audit's contig-length assertions change; the
translation logic is species-neutral.

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

## Checking the reference itself

```bash
v2p audit --ref ref/
```

Verifies that the GTF, genome and proteome are being *read* correctly,
separately from whether the output looks right. Four real bugs were found
this way that output-agreement testing missed entirely. Run it whenever you
change reference release.

# Conversion plan — HCC1395 variants → amino-acid FASTA

## 0. Target format (now confirmed)

`Human_Homo_sapiens_uniprot_SP_UP000005640_N20432_20260304.fasta` is the
spec Zhiyin meant. Measured properties:

| property | value |
|---|---|
| entries | 20,432 |
| database | SwissProt only (`sp\|`), no TrEMBL |
| isoforms | none — canonical only, no `-2` accessions |
| header | `>sp\|ACC\|ENTRY_HUMAN Description OS=Homo sapiens OX=9606 [GN=SYM]` |
| fields present | `OS=` 20,432 · `OX=` 20,432 · `GN=` 20,307 (125 lack `GN=`) |
| fields absent | `PE=`, `SV=` — do **not** add them |
| sequence layout | **one line per entry, unwrapped** (all 20,432) |
| alphabet | 20 standard + `U` (36 selenocysteine residues) |
| length | min 2, median 415, max 34,350 |

The `uniprot` output style reproduces this exactly. Variant entries use a
`vr|` prefix so a search engine can separate variant from canonical hits
on the database field alone; variant metadata is appended as further
`KEY=value` pairs, which is how UniProt already extends the line.

```
>vr|HCC1395_RNAEDIT_000007|GRIA2_HUMAN_RNA_EDITING Glutamate receptor 2 (p.Q607R; RNA editing site (A-to-I)) OS=Homo sapiens OX=9606 GN=GRIA2 VT=RNA_EDITING CSQ=missense TX=NM_000826 PC=p.Q607R LOC=chr4:157336723A>G POS=607 NOVEL=607-607 SRC=RES_ANNOVAR QC=ref_residue_match
MQKIMHISVLLSPVLWGLIFGVSSNSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTP...
```

---

## 1. Workflow

```
        ┌─ VCF (sSNV/sIndel) ──┐
        ├─ RES multianno ──────┤
INPUTS  ├─ Fusion CSV ─────────┤──► [1] parse ──► unified manifest (9,148 + VCF)
        └─ AS-LR CSV ──────────┘                        │
                                                        ▼
REFERENCE  GRCh38 + GENCODE v44 ──────────────► [2] translate
           UniProt SwissProt ──┐                   ├── smallvar  (SNV/InDel/RES)
                               │                   ├── fusion    (chimeric tx)
                               │                   └── splicing  (isoform)
                               │                        │
                               │                        ▼
                               ├──────────────► [4] UniProt validation gate
                               │                   A. annotation residue check
                               │                   B. reference-translation check
                               │                        │
        ORTHOGONAL TOOLS       │                        ▼
        VEP ProteinSeqs ───────┼──────────────► [5] tool comparison
        AGFusion ──────────────┤                   sequence + 9-mer Jaccard
        gffread -y ────────────┘                        │
                                                        ▼
                                              [6] combined UniProt-style FASTA
                                                  + record table + logs
```

Stages 1, 3 and 4A are reference-free and **have already been run** on your
data. Stages 2, 4B, 5 and 6 need the reference downloads.

---

## 2. Quality control — seven gates

Not opinions about the code; each one is a check that can fail and stop the
run.

| # | gate | what it catches | status |
|---|---|---|---|
| **G1** | Unit + integration tests on a synthetic genome | translation, strand, frame, indel and fusion logic errors | **29/29 pass** |
| **G2** | Input consistency audit | corrupt or mislabelled inputs — e.g. all 8,093 RES sites are A>G/T>C, exactly as A-to-I editing predicts | **run, clean** |
| **G3** | REF-allele verification per variant | wrong genome build, wrong coordinate convention, off-by-one | in stage 2, flags `REF_MISMATCH` |
| **G4a** | UniProt annotation check | isoform/transcript mismatch between the caller and the reference proteome | **run: 19/25 sites anchored** |
| **G4b** | UniProt translation check | our reference proteins vs UniProt residue-by-residue; **fails the run below 90% agreement** | pending stage 2 |
| **G5** | Orthogonal tool comparison | systematic errors a single implementation can't see | script ready |
| **G6** | Biological positive controls | known-answer sanity check | **passed** — GRIA2 Q607R, NEIL1 K242R, BLCAP Y2C, CDK13 Q103R, COG3 I635V all recovered, all canonical ADAR recoding sites |
| **G7** | Provenance | irreproducibility | SHA-256 + `pip freeze` + full params per stage |

**Why G4b is the strongest gate.** If our translated reference GRIA2 does
not equal UniProt's GRIA2, the variant GRIA2 is worthless. This validates
the genome build, the GENCODE release, the exon assembly, the CDS offset,
the strand handling and the codon table in one comparison, against a
reference you supplied — no second tool required. Selenoproteins (25 genes)
are excluded from the statistic rather than counted as failures, since a
standard-code translation legitimately truncates at UGA.

**G4a result already in hand.** Of 25 unique recoding sites: 19 anchored,
2 genes absent from SwissProt (RPSAP58 is a pseudogene, TMEM183B a
paralog), 4 with no agreeing transcript (AASDH, CCDC144A, SRP9, TUBGCP2 —
each annotated only on a non-canonical RefSeq isoform). The pattern is
consistent: wherever ANNOVAR reports several transcripts, at least one
matches UniProt canonical numbering. That is isoform choice, not a wrong
call — but every one is listed in `tables/uniprot_residue_check.tsv`.

---

## 3. Tools

Our pipeline is primary because no single external tool covers all four
variant classes. Each class then gets an independent second opinion.

| variant class | primary | orthogonal cross-check | why that tool |
|---|---|---|---|
| SNV / InDel | `v2p.build.smallvar` | **Ensembl VEP + `ProteinSeqs` plugin** | the reference implementation; emits reference and mutated protein FASTA directly |
| SNV / InDel (2nd) | — | **pypgatk `vcf-to-proteindb`** | pure-Python, Ensembl-based, easy to diff |
| RNA editing | `v2p.build.smallvar` | **VEP**, after converting the RES table to VCF | an A-to-I edit is a transcript-level A>G; VEP handles it if you hand it a VCF |
| Fusion | `v2p.build.fusion` | **AGFusion** | purpose-built; `--middlestar` marks the junction in the peptide |
| Splicing | `v2p.build.splicing` | **`gffread -y`** on our emitted isoform GTF | independent CDS→protein translator; isolates our translation from our isoform construction |
| all | — | **SnpEff** | second opinion on consequence class only |

Deliberately **not** used: `customProDB` (unmaintained, R-only, Ensembl-89
era annotation) and full neoantigen pipelines like pVACtools/nextNEOpi
(they predict MHC binding, which is a different question from "convert
these variants to amino-acid sequences" and would add HLA typing you don't
need).

---

## 4. Install

```bash
# --- core (needed for our pipeline) --------------------------------------
conda create -n hcc1395 -c conda-forge -c bioconda python=3.11 \
    pyfaidx biopython samtools bcftools gffread
conda activate hcc1395

# --- orthogonal tools ----------------------------------------------------
conda install -c bioconda ensembl-vep=112 snpeff      # VEP + SnpEff
pip install agfusion pypgatk                          # fusions + 2nd SNV tool

# --- VEP cache and the ProteinSeqs plugin --------------------------------
vep_install -a cfp -s homo_sapiens -y GRCh38 -c $HOME/.vep \
            -g ProteinSeqs,Wildtype,Frameshift

# --- AGFusion database (must match the GENCODE/Ensembl release you use) ---
agfusion download -g hg38 -r 112
```

Rough disk: VEP cache ~26 GB, GRCh38 FASTA ~3 GB, AGFusion DB ~2 GB.

If conda is awkward, VEP also runs from Docker:
`docker run -v $PWD:/data ensemblorg/ensembl-vep:release_112 vep ...`

---

## 5. Commands

### 5.0 Inspect everything first

```bash
unzip -o High_conf_variant_dataset.zip -d data/
python3 scripts/00_inspect_inputs.py data/ --json
```

Identifies each file by content rather than extension and reports what is
actually inside: for a VCF the sample columns, FILTER tiers, variant-type
breakdown and **the genome build read from the contig lengths**; for a BED
the region count and genome fraction; for a FASTA the header grammar,
wrapping and alphabet. It then cross-checks contig naming and build across
all files and refuses to be quiet about a mismatch. Standard library only,
so it runs before you install anything.

The `.csi` needs no action — it is an index for the bgzipped VCF and can be
regenerated at any time with `bcftools index`.

### 5.1 References

```bash
bash scripts/00_fetch_references.sh ref/
cp <your file> ref/uniprot_human_SP.fasta
```
The script asserts `chr1 == 248,956,422 bp` and aborts if you have GRCh37.

### 5.2 Our pipeline

```bash
python scripts/01_parse_inputs.py \
    --vcf    data/high-confidence_sSNV_sIndel_v1.sort.final.combined.sort.vcf \
    --res    data/HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt \
    --fusion data/HCC1395_high_confidence_Fusion_genes_all.csv \
    --as-lr  data/HCC1395_high_confidence_AS-LR_v1.csv \
    --outdir results

python scripts/03_qc_report.py --manifest results/tables/unified_variant_manifest.tsv

python scripts/02_build_protein_fasta.py \
    --manifest results/tables/unified_variant_manifest.tsv \
    --genome   ref/GRCh38.primary_assembly.genome.fa \
    --gtf      ref/gencode.v44.annotation.gtf.gz \
    --uniprot  ref/uniprot_human_SP.fasta \
    --header-style uniprot --include-reference \
    --outdir   results

python scripts/04_validate_uniprot.py \
    --uniprot ref/uniprot_human_SP.fasta \
    --recoding results/tables/res_recoding_sites.tsv \
    --protein-fasta results/fasta/HCC1395_variant_proteins.uniprot.fasta
```

Stage 4 exits non-zero if reference-translation agreement falls below 90%.
Wire that into CI and a bad reference build can never ship.

### 5.3 Cross-check A — VEP on SNVs/InDels

```bash
vep --input_file data/high-confidence_sSNV_sIndel_v1.sort.final.combined.sort.vcf \
    --output_file vep_out/annotated.txt \
    --species homo_sapiens --assembly GRCh38 --cache --offline \
    --dir_cache $HOME/.vep --fasta ref/GRCh38.primary_assembly.genome.fa \
    --symbol --protein --uniprot --hgvs --coding_only --pick \
    --plugin ProteinSeqs,vep_out/reference.fa,vep_out/mutated.fa
```

`ProteinSeqs` disables `--fork`, so this is single-threaded and slow on
~40k variants — budget an hour. `--pick` gives one transcript per variant,
matching our default `--transcript-mode representative`; drop it and use
`--transcript-mode all` on our side if you want the comparison to be
like-for-like across every transcript.

### 5.4 Cross-check B — VEP on RNA editing

VEP needs a VCF, so convert the ANNOVAR table first:

```bash
python - <<'PY'
import csv
rows = csv.DictReader(open('data/HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt'), delimiter='\t')
with open('res_sites.vcf','w') as o:
    o.write('##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for r in rows:
        if not r['Chr'] or r['Chr'].upper()=='CHROM': continue
        o.write(f"{r['Chr']}\t{r['Start']}\t.\t{r['Ref']}\t{r['Alt']}\t.\tPASS\tSOURCE=RES\n")
PY
bcftools sort res_sites.vcf -o res_sites.sorted.vcf
vep -i res_sites.sorted.vcf -o vep_out/res.txt --cache --offline --assembly GRCh38 \
    --fasta ref/GRCh38.primary_assembly.genome.fa --symbol --protein --hgvs --pick \
    --plugin ProteinSeqs,vep_out/res_reference.fa,vep_out/res_mutated.fa
```

Expected: VEP reports ~25 nonsynonymous sites, matching ANNOVAR. A large
discrepancy means the two annotation sources disagree on transcript models,
which you need to know before shipping.

### 5.5 Cross-check C — AGFusion on fusions

```bash
mkdir -p agf_out
tail -n +2 data/HCC1395_high_confidence_Fusion_genes_all.csv | tr -d '\r' | \
while IFS=, read -r tag g1 g2 type validated conf bp1 bp2; do
  agfusion annotate \
      --gene5prime "$g1" --junction5prime "${bp1##*:}" \
      --gene3prime "$g2" --junction3prime "${bp2##*:}" \
      -db agfusion.homo_sapiens.112.db -o "agf_out/${g1}--${g2}" \
      --middlestar --noncanonical 2>>agf_out/errors.log || true
done
cat agf_out/*/*_protein.fa > agf_out/all_fusion_proteins.fa
```

`--middlestar` inserts a `*` at the junction, so you can confirm AGFusion
and we place the junction at the same residue. Expect failures on the
lncRNA partners (`RP11-...`, `CTD-...`, `AC0...`) — AGFusion needs
protein-coding genes, and 8 of your 28 fusions have a non-coding partner.
Those are legitimately out of scope for both tools.

### 5.6 Cross-check D — gffread on splicing isoforms

```bash
python scripts/02_build_protein_fasta.py \
    --manifest results/tables/unified_variant_manifest.tsv \
    --genome ref/GRCh38.primary_assembly.genome.fa \
    --gtf ref/gencode.v44.annotation.gtf.gz \
    --classes AS_SE,AS_RI,AS_A3,AS_A5,AS_MX,AS_AF,AS_AL \
    --emit-isoform-gtf results/tables/as_isoforms.gtf \
    --outdir results/as_only

gffread -y gffread_out/as_proteins.fa \
        -g ref/GRCh38.primary_assembly.genome.fa \
        results/tables/as_isoforms.gtf
```

This separates the two things that can go wrong in the splicing path: if
gffread and we agree, our translation is right and any remaining doubt is
about isoform *construction*, which is a modelling choice, not a bug.

### 5.7 Compare everything

```bash
python scripts/05_compare_tools.py \
    --fasta ours=results/fasta/HCC1395_variant_proteins.uniprot.fasta \
    --fasta vep=vep_out/mutated.fa \
    --fasta agfusion=agf_out/all_fusion_proteins.fa \
    --fasta gffread=gffread_out/as_proteins.fa \
    --k 9 --outdir results
```

Reports, per pair: shared complete sequences, Jaccard on sequences, shared
9-mers, Jaccard on 9-mer space, and the list of genes both tools
translated but on which they disagree.

**How to read it.** High k-mer Jaccard with low sequence Jaccard means the
tools agree on the peptides and differ only on transcript choice — decide a
policy and move on. Low k-mer Jaccard is a real disagreement and needs
gene-by-gene inspection. 9-mer space is the right unit because it is what
determines whether a downstream MS search or binding predictor can find a
peptide at all.

**Acceptance thresholds** (proposed, adjust with Zhiyin):

| comparison | expected 9-mer Jaccard | action if below |
|---|---|---|
| ours vs VEP (SNV/InDel, `--pick` both) | > 0.95 | inspect disagreeing genes individually |
| ours vs VEP (RNA editing) | > 0.95 | check ANNOVAR vs Ensembl transcript models |
| ours vs AGFusion (coding partners only) | > 0.80 | compare junction residue placement |
| ours vs gffread (splicing) | > 0.99 | our translation is wrong, not the isoform model |

---

## 6. Still outstanding

1. **The SNV/InDel VCF.** The largest class (~39k SNVs, ~1.6k indels) and
   the only one not yet received. Everything else is in place and waiting
   for it.
2. **`PGx_High-Confidence_Regions_v1.6.sorted.bed`.** Optional but useful:
   restricting to high-confidence regions cuts false positives. Add
   `--regions <bed>` support once you send it.
3. **Confirm with Zhiyin:** should the delivered file be variant entries
   only, or variant entries **concatenated with** the 20,432 canonical
   SwissProt sequences? For an MS database search you want both — a
   variant peptide is only interpretable next to its wild-type counterpart.
   `--include-reference` currently emits only the reference proteins of
   *affected* transcripts, not the whole proteome; say the word and it can
   append the full 20,432.
4. **Transcript policy — resolved: build both.** `make both` runs the same
   manifest through `--transcript-mode representative` and
   `--transcript-mode all`, then compares the two databases with stage 5.
   Only that one flag differs, so the reported 9-mer Jaccard *is* the cost
   of the choice. If the two cover nearly the same peptide space, ship the
   representative database — it is 4–6× smaller and gives a better FDR on
   any downstream search. If they diverge, the extra transcripts are
   contributing real sequence and are worth their size.

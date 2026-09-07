# Output FASTA format specification

> **Open item.** Zhiyin's email refers to a target FASTA format
> specification "attached for reference", but no specification file arrived
> with the package (the three data files plus the SNV/InDel screenshot were
> all that came through). Rather than guess, the writer supports three
> conventions behind a `--header-style` flag. Switching is a one-flag
> re-emission and does **not** require re-running translation, so if the
> real specification turns out to be a fourth layout, adding it is a single
> formatter function in `src/v2p/fasta.py` — the sequences are unaffected.

All three styles carry the variant-type information the email asked for.

---

## Variant-type vocabulary

Emitted verbatim in every style. Stable across releases.

| value | meaning |
|---|---|
| `SNV` | somatic single-nucleotide variant |
| `MNV` | somatic multi-nucleotide variant |
| `INDEL` | somatic insertion / deletion |
| `RNA_EDITING` | RNA editing site (A-to-I) |
| `FUSION` | gene fusion |
| `AS_SE` | alternative splicing: skipped exon |
| `AS_RI` | alternative splicing: retained intron |
| `AS_A3` | alternative splicing: alternative 3′ splice site |
| `AS_A5` | alternative splicing: alternative 5′ splice site |
| `AS_MX` | alternative splicing: mutually exclusive exons |
| `AS_AF` | alternative splicing: alternative first exon |
| `AS_AL` | alternative splicing: alternative last exon |
| `REFERENCE` | unmodified reference protein (with `--include-reference`) |

## Consequence vocabulary

`synonymous`, `missense`, `stop_gained`, `stop_lost`, `start_lost`,
`inframe_insertion`, `inframe_deletion`, `frameshift`, `in_frame_fusion`,
`frameshift_fusion`, `5_prime_UTR`, `3_prime_UTR`,
`<TYPE>_form1_annotated_isoform`, `<TYPE>_form2_annotated_isoform`,
`<TYPE>_constructed_isoform`, `reference`.

`5_prime_UTR` and `3_prime_UTR` are retained under `--keep-unchanged` (the
default): the variant is real and lies in the transcript, but outside the
CDS, so the protein is identical to the wild type. They carry their own
labels rather than being folded into `synonymous`, because conflating the
two once inverted the nonsynonymous:synonymous ratio and shipped sixteen
indels as false frameshifts.

---

## Style 1 — `peff` (recommended)

PSI Extended FASTA Format, the ratified HUPO-PSI standard for encoding
sequence variants in FASTA. Backward compatible: a plain FASTA parser reads
it as an ordinary file and ignores the keys. Natively understood by Comet,
neXtProt and UniProtKB.

File header block:

```
# PEFF 1.0
# //
# DbName=HCC1395_variant_proteome
# DbDescription=HCC1395 variant protein database
# Prefix=HCC1395
# DbSource=HCC1395 high-confidence variant package
# DbVersion=<YYYY-MM-DD>
# SequenceType=AA
# NumberOfEntries=<n>
# GeneralComment=Custom keys ... are non-standard extensions; PEFF-compliant
#                readers ignore unknown keys.
# //
```

Entry line keys:

| key | standard? | content |
|---|---|---|
| `\DbUniqueId` | PSI | unique entry id |
| `\PName` | PSI | gene + variant-type description |
| `\GName` | PSI | gene symbol |
| `\TaxName` / `\NcbiTaxId` | PSI | `Homo sapiens` / `9606` |
| `\Length` | PSI | residue count |
| `\VariantSimple` | PSI | `(pos\|newAA)` — emitted for missense and stop-gain only |
| `\VariantType` | **custom** | variant-type vocabulary above |
| `\Consequence` | **custom** | consequence vocabulary above |
| `\TranscriptId` | **custom** | transcript(s) used |
| `\ProteinChange` | **custom** | HGVS-like, e.g. `p.Q607R`, `p.L6Ffs*17`, `p.junction@11` |
| `\GenomicLocus` | **custom** | e.g. `chr4:157336723A>G` or `chr5:78267520::chr5:74768931` |
| `\NovelSpan` | **custom** | 1-based inclusive residue range that differs from reference |
| `\EvidenceSource` | **custom** | which input file the record came from |
| `\Confidence` | **custom** | source confidence (VCF FILTER, fusion score, AS confidence) |
| `\Comment` | PSI | semicolon-joined processing notes / warnings |

`\VariantSimple` has no representation for frameshifts, fusion junctions or
whole isoforms — those are complete alternative sequences, not residue
substitutions — so for those classes the type is carried by `\VariantType`
and the entry is a standalone sequence.

## Style 2 — `descriptive`

Space-separated `KEY=value` after a pipe-delimited id. Safest with search
engines that mangle backslashes.

```
>HCC1395|<seq_id> VT=<type> CSQ=<consequence> GN=<gene> TX=<transcript>
 PC=<protein change> LOC=<locus> SRC=<source> CONF=<confidence> LEN=<n>
 [POS=<aa pos>] [NOVEL=<start>-<end>] [NOTE=<comma-joined notes>]
```

## Style 3 — `pvac`

Compact, dot-delimited, for neoantigen tooling that expects `MT.`/`WT.` prefixes:

```
>MT.<gene>.<transcript>.<variant type>.<consequence>[.<change>]
>WT.<gene>.<transcript>.REFERENCE.reference
```

---

## Companion table

`results/tables/protein_records.tsv`, one row per FASTA entry:

`seq_id`, `variant_class`, `consequence`, `gene`, `transcript`,
`protein_change`, `variant_pos_aa`, `protein_length`, `ref_protein_len`,
`locus`, `source`, `confidence`, `novel_span`, `notes`

## Processing notes vocabulary

Flags that appear in `\Comment` / `NOTE=` / the `notes` column:

| note | meaning |
|---|---|
| `REF_MISMATCH(expected=X,found=Y)` | the call's REF allele disagrees with the reference — **check the build** |
| `5UTR_indel_shifts_CDS_start` | an indel upstream of the start codon moved the CDS |
| `alt_orf_runs_to_transcript_end` | no stop codon before the transcript ends |
| `ambiguous_bases_translated_as_X` | reference `N` run inside the ORF |
| `junction_on_codon_boundary` | fusion: 5′ partner frame intact at the junction |
| `3p_native_frame` / `3p_frameshifted` | fusion: 3′ partner in its own frame or not |
| `bp5_intronic_snapped_to_exon_end` | fusion breakpoint moved to the splice donor |
| `bp3_intronic_snapped_to_exon_start` | fusion breakpoint moved to the splice acceptor |
| `matched_annotated_transcript` | AS form found in GENCODE — high confidence |
| `constructed_isoform` | AS form built by editing the exon chain — lower confidence |
| `start_codon_inferred_first_ATG` | annotated start unavailable; first viable ATG used |

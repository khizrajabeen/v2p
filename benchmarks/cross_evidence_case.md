# A cross-evidence proteoform, constructed

**This is a constructed case, not a measurement.** One half is a real call from HCC1395; the other is synthetic, and which is which is stated below.

## Why it had to be constructed

HCC1395 contains no transcript where an RNA editing site and a DNA variant both recode. That was measured, not assumed:

| check | count |
|---|---|
| transcripts carrying both an editing site and a DNA variant | 40 |
| of those, transcripts where **both** recode | **0** |
| genes with a recoding edit *and* a recoding DNA variant | **0** |

29 genes carry a recoding edit and 278 carry a recoding DNA variant; the two sets do not intersect.

## The case

Gene **POTEE**, transcript `ENST00000683005.1`.

| element | status | detail |
|---|---|---|
| A-to-I editing site | **real**, called in HCC1395 | `chr2:131250822A>G`, residue 509 |
| somatic SNV | **synthetic** | `chr2:131245453T>G`, residue 486 |

The editing site is the half no genotype-based tool can use. The SNV stands in for a somatic call this sample does not have at that position.

## What the entries contain

| entry | residue 486 | residue 509 |
|---|---|---|
| reference | Y | K |
| edit only | Y | E |
| SNV only | D | K |
| **combined** | D | E |

v2p describes the combined entry as `p.[Y486D(;)K509E]`, class `COMBO`, notes `combinatorial,unphased,cross_evidence_RNA_EDITING+SNV`.

## The peptide that exists only in the combined form

```
DSSENSNPEQDLKLTSEEESQRLEGSENGQPEK
```

1 tryptic peptide(s) of the combined protein appear in neither single-variant entry nor the reference. A search against a database built one variant at a time cannot identify them at any FDR.

## What this does and does not show

It shows the mechanism, and exercises the real code path on a real transcript with a real editing site. It does **not** show how often such proteoforms occur in nature. That needs a sample with both call types co-located, which the search recorded in `benchmarks/README.md` did not find.

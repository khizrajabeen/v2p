# Example input folder

A small input folder so `v2p detect` and `v2p run` work straight after a
clone. Point the tool at this directory:

```bash
v2p detect examples/
v2p run examples/ --ref ref/ --outdir out/ --name EXAMPLE --logdir logs/
```

| file | evidence type | rows |
|---|---|---|
| `small_variants.vcf` | somatic SNV / InDel | 8 |
| `rna_editing.txt` | A-to-I RNA editing (ANNOVAR-style) | 5 |
| `fusions.csv` | gene fusion calls | 3 |
| `splicing.csv` | alternative splicing (SUPPA2 event ids) | 5 |

## What these are, and are not

These are **synthetic calls at real coordinates**. Nothing here is
redistributed from the HCC1395 high-confidence package, which is not ours
to publish.

- **Small variants** are placed in the CDS of eight well-known genes
  (TP53, KRAS, BRAF, PIK3CA, EGFR, PTEN, NRAS, IDH1), at a different codon
  in each. The reference allele is read back from the genome, so each call
  is internally consistent and produces a genuine missense consequence.
  They are *not* real observed mutations — do not treat them as such.
- **RNA editing** rows are the canonical ADAR recoding sites: GRIA2 Q607R,
  NEIL1 K242R, BLCAP Y2C, CDK13 Q103R, COG3 I635V. These are real,
  published recoding events, and they double as the pipeline's positive
  controls: a run over this folder should recover each one.
- **Fusions** are published, characterised rearrangements — BCR–ABL1,
  EML4–ALK, TMPRSS2–ERG — at their canonical breakpoints.
- **Splicing** events are built from real GENCODE v44 junctions.

## Regenerating

```bash
python examples/make_examples.py --genome ref/... --gtf ref/...
```

Needs the references, so it is not part of `make test`. Re-run it only if
the coordinates or the annotation release change.

The ADAR coordinates are **derived, not quoted**. For each published
protein change the generator locates that codon in the GENCODE
representative transcript, and writes the site only if a single A→G in
that codon reproduces the published substitution *and* the genome base
matches the expected strand. A site that cannot be derived is reported and
dropped rather than shipped as a guess.

That check earns its keep: four of the five ADAR coordinates originally
written from memory were wrong, and the derivation caught all four.

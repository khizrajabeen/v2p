| measure | v2p | ProVar |
|---|---|---|
| variants represented in >=1 entry | 2865 | 911 |
| sequences produced | 4260 | 2197 |
| entries carrying >1 variant | 7 | 0 |
| novel tryptic peptides | 207411 | 310372 |
| peptides the other tool lacks | 133229 | 236190 |
| peptides from combined entries | 121 | 0 |
| wall clock | 2:18.50 | 23:55.70 |
| peak memory | 472 MB | 839 MB |

Variants ProVar represents and v2p does not: **8**

Each of these is a v2p bug until investigated. What the investigation found is in `README.md`.


- `11:61740459TCAGCTGCTGTGAGCAGGAGGAGCCCACATACTTTGC>T`
- `12:104251652AGGTCCACATGCACACGCTGTACTGAGGTAAGGCTTTAAACTCAAGG>A`
- `13:36972022TATAGCCATTTGATTCTGAGGGAAAAAGCAGAG>T`
- `19:37692270G>A`
- `1:155086199GTGCTGGGTGAGTCTGCGCAGCGCCCTCTGGTGGCCAC>G`
- `1:235766274ATTCCTGAAAAAATAAAAAAAACTCTCT>A`
- `2:26182470CGATGCACAGGATGCCAGATCCAGGTAGGG>C`
- `3:194137671TCTGGCCCGCAGCTGCGCTGAGCACAGACCCAA>T`

Variants v2p represents and ProVar does not: **1962**

## Caveats, which the numbers above do not carry on their own

- **The peptide counts are confounded.** Novelty here is measured against UniProt SwissProt, the reference v2p appends. ProVar translates Ensembl transcripts, so its alternative isoforms count as novel whether or not a variant is involved. That inflates both its novel-peptide total and its peptides-the-other-lacks column, and it is why ProVar can show more novel peptides from fewer sequences. Read those two rows as *different reference sets*, not as variant discovery.
- **Variant coverage is the clean comparison**, because a variant either ends up in an entry or it does not, independent of which reference proteome either tool started from.
- **The runtime figures are not like for like.** This harness re-reads the 521 MB cDNA FASTA once per chromosome, which is most of ProVar's wall clock. It is a property of how it was driven here, not a claim about the tool.
- **Entries carrying more than one variant is the row that matters** for the feature under test, and it is measured, not assumed: ProVar's own output confirms one variant per entry.

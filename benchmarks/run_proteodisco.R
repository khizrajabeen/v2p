#!/usr/bin/env Rscript
# Run ProteoDisco over the same truth VCF that v2p and pypgatk were given.
#
# The published recipe uses BSgenome.Hsapiens.UCSC.hg38 and
# TxDb.Hsapiens.UCSC.hg38.knownGene. This deliberately does not. Those are
# UCSC knownGene, a different transcript set from the GENCODE v44 v2p
# uses, so a comparison built on them measures two things at once: the
# tools, and the annotations they were handed. Building the TxDb from our
# own GTF and reading the genome from our own FASTA gives every tool the
# same reference, which is the only way the numbers mean anything.
#
# It is also far lighter: those two annotation packages are together over
# a gigabyte, and the GTF and genome are already on disk.
#
#   Rscript benchmarks/run_proteodisco.R \
#       <truth.vcf> <gencode.gtf.gz> <genome.fa> <out.fasta>

suppressPackageStartupMessages({
  library(ProteoDisco)
  library(GenomicFeatures)
  library(Rsamtools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: run_proteodisco.R <vcf> <gtf.gz> <genome.fa> <out.fasta>")
}
vcf_path <- args[1]
gtf_path <- args[2]
fa_path  <- args[3]
out_path <- args[4]

message("building a TxDb from ", gtf_path, " (several minutes) ...")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")

message("opening the genome ...")
genome <- Rsamtools::FaFile(fa_path)

message("generating the ProteoDiscography ...")
pd <- ProteoDisco::generateProteoDiscography(
  TxDb = txdb,
  genomeSeqs = genome
)

message("importing variants ...")
pd <- ProteoDisco::importGenomicVariants(
  pd,
  files = vcf_path,
  samplenames = "truth"
)

message("incorporating ...")
pd <- ProteoDisco::incorporateGenomicVariants(
  pd,
  aggregateSamples = FALSE,
  aggregateWithinTranscript = FALSE
)

message("writing ", out_path, " ...")
ProteoDisco::exportProteoDiscography(pd, outFile = out_path)

# The comparison scores locus coverage; the Python side does the scoring
# against the same truth set, so just report what was produced.
res <- ProteoDisco::mutantTranscripts(pd)
n <- if (is.null(res$genomicVariants)) 0 else nrow(res$genomicVariants)
message("done. variant transcripts: ", n)

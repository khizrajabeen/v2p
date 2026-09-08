# HCC1395 variant -> protein FASTA
#
#   make test                                  run the test suite (no reference needed)
#   make parse                                 stage 1 + 3 (no reference needed)
#   make all GENOME=... GTF=... [VCF=...]      full pipeline
#
GENOME ?= ref/GRCh38.primary_assembly.genome.fa
GTF    ?= ref/gencode.v44.annotation.gtf.gz
VCF    ?=
STYLE  ?= peff
OUT    ?= results
DATA   ?= data
MANIFEST = $(OUT)/tables/unified_variant_manifest.tsv
VCFARG = $(if $(VCF),--vcf $(VCF),)

.PHONY: all test reproducibility refs parse fasta qc clean
all: parse fasta

test:
	python3 tests/test_pipeline.py
	python3 tests/test_invariants.py
	python3 tests/test_config.py
	python3 tests/test_species.py
	python3 tests/test_noncanonical.py
	python3 tests/test_combinatorial.py

# M2 acceptance test. Needs the full reference, so it is deliberately not
# part of `test`: two runs from one config must produce a byte-identical
# release, which is what catches set-iteration order, unsorted globs and
# timestamps or absolute paths embedded in output files.
REPRO_DIR ?= .repro
reproducibility:
	@rm -rf $(REPRO_DIR)
	python3 scripts/v2p run data --ref ref --name REPRO \
	    --write-config $(REPRO_DIR)/run.yaml
	python3 scripts/v2p run --config $(REPRO_DIR)/run.yaml \
	    --outdir $(REPRO_DIR)/a --logdir $(REPRO_DIR)/log_a
	python3 scripts/v2p run --config $(REPRO_DIR)/run.yaml \
	    --outdir $(REPRO_DIR)/b --logdir $(REPRO_DIR)/log_b
	@diff -r $(REPRO_DIR)/a $(REPRO_DIR)/b \
	  && echo "PASS: two runs from one config are byte-identical" \
	  || (echo "FAIL: the two runs differ"; exit 1)

refs:
	bash scripts/00_fetch_references.sh ref/

parse: $(MANIFEST) qc

$(MANIFEST):
	python3 src/v2p/stages/01_parse_inputs.py $(VCFARG) \
	  --res    $(DATA)/HCC1395_high_confidence_RES_v1_addAlu_hg38_multianno.txt \
	  --fusion $(DATA)/HCC1395_high_confidence_Fusion_genes_all.csv \
	  --as-lr  $(DATA)/HCC1395_high_confidence_AS-LR_v1.csv \
	  --outdir $(OUT)

qc: $(MANIFEST)
	python3 src/v2p/stages/03_qc_report.py --manifest $(MANIFEST) --outdir $(OUT)

MODE  ?= representative

fasta: $(MANIFEST)
	python3 src/v2p/stages/02_build_protein_fasta.py \
	  --manifest $(MANIFEST) --genome $(GENOME) --gtf $(GTF) \
	  --uniprot $(UNIPROT) --header-style $(STYLE) --transcript-mode $(MODE) \
	  --include-reference --outdir $(OUT)

# Build both transcript policies and measure how much the choice costs you.
# Same manifest, same reference, only --transcript-mode differs, so any
# difference in the comparison is attributable to that flag alone.
both: $(MANIFEST)
	$(MAKE) fasta MODE=representative
	$(MAKE) fasta MODE=all
	python3 src/v2p/stages/05_compare_tools.py \
	  --fasta representative=$(OUT)/fasta/HCC1395_variant_proteins.$(STYLE).representative.fasta \
	  --fasta all_transcripts=$(OUT)/fasta/HCC1395_variant_proteins.$(STYLE).all.fasta \
	  --k 9 --outdir $(OUT)

clean:
	rm -rf $(OUT)/fasta $(OUT)/tables $(OUT)/qc $(OUT)/_itest logs/*.log logs/*.json

# ---- validation and cross-checks ----------------------------------------
UNIPROT ?= ref/uniprot_human_SP.fasta

.PHONY: validate compare
validate:
	python3 src/v2p/stages/04_validate_uniprot.py --uniprot $(UNIPROT) \
	  --recoding $(OUT)/tables/res_recoding_sites.tsv \
	  $(if $(wildcard $(OUT)/fasta/*.fasta),--protein-fasta $(firstword $(wildcard $(OUT)/fasta/*.fasta)),) \
	  --outdir $(OUT)

compare:
	python3 src/v2p/stages/05_compare_tools.py \
	  --fasta ours=$(OUT)/fasta/HCC1395_variant_proteins.uniprot.fasta \
	  --fasta vep=vep_out/mutated.fa \
	  --fasta agfusion=agf_out/all_fusion_proteins.fa \
	  --k 9 --outdir $(OUT)

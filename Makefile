
.PHONY: all
all: report

nextflow_path = scripts/nextflow

tmp/counts_table.csv: $(nextflow_path) \
          scripts/analyze_barcode_sequencing.nf \
          scripts/analyze_barcode_sequencing.nfconfig
	$< run scripts/analyze_barcode_sequencing.nf -C analyze_barcode_sequencing.nfconfig


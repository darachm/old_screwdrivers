
.PHONY: all
all: report

nextflow_path = scripts/nextflow

tmp/counts_table.csv: $(nextflow_path) \
          scripts/quantify_barcodes.nf \
          scripts/quantify_barcodes.nfconfig
	$< run scripts/quantify_barcodes.nf -C quantify_barcodes.nfconfig


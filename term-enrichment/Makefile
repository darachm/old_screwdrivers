
.PHONY: all
all: term_enrichment_clusterProfiler.html

data_files = data/Figure4_Table_BFFfilteredPooledModels.csv \
  data/sgd_go_slim_171013.txt data/sgd_go_full_171013.txt \
  data/sgd_go_terms_171013.txt

term_enrichment_clusterProfiler.html \
  output_table_categorical.csv output_table_GSEA.csv \
  : \
  term_enrichment_clusterProfiler.Rmd $(data_files)
	Rscript -e "rmarkdown::render('$<')"


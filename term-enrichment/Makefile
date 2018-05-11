
.PHONY: all
all: term_enrichment_clusterProfiler.html

term_enrichment_clusterProfiler.html \
  : \
  term_enrichment_clusterProfiler.Rmd 
	Rscript -e "rmarkdown::render('$<')"



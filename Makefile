.PHONY: all
all: something_primitive.pdf

something_primitive.pdf: something_primitive.Rmd preamble-latex.tex
	Rscript --vanilla -e "rmarkdown::render('something_primitive.Rmd')"



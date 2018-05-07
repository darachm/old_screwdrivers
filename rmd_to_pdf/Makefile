.PHONY: all
all: something.pdf something.html


something.pdf something.html \
  : \
  something_primitive.Rmd scripts/preamble-latex.tex
	Rscript --vanilla -e "rmarkdown::render('$<')"

#,output_file='../output/Figure2_Supplementary_Writeup.pdf',output_dir='output/')"

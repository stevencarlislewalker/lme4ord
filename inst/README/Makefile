all: README.md copy

README.md: README.Rmd
	Rscript -e "library(knitr); knit(\"README.Rmd\")"
	rm -r -f figure/
	mv inst/README/figure/ figure/
	rm -r inst

copy: README.md
	cp README.md ../../README.md
	rm README.md

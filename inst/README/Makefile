all: README.md copy

README.md: README.Rmd
	Rscript -e "library(knitr); knit(\"README.Rmd\")"
	make clean
	mv inst/README/figure/ figure/
	rm -r inst
	git add figure ../../README.md README.Rmd

copy: README.md
	cp README.md ../../README.md
	rm README.md

README.html: README.md
	pandoc README.md -o README.html

README.R: README.Rmd
	Rscript -e "library(knitr); purl(\"README.Rmd\")"

clean:
	git rm -r -f figure/

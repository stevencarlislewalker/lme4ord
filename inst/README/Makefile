all: README.md copy git

README.md: README.Rmd
	Rscript -e "library(knitr); knit(\"README.Rmd\")"
	make rmfig
	mv inst/README/figure/ figure/
	rm -r inst

copy: README.md
	cp README.md ../../README.md
	rm README.md

README.html: README.md
	pandoc README.md -o README.html

README.R: README.Rmd
	Rscript -e "library(knitr); purl(\"README.Rmd\")"

git: 
	git add figure ../../README.md README.Rmd

rmfig:
	git rm -r -f figure

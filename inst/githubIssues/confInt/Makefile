SRC      := confInt
PKGS     := --default-packages=base,methods,stats,utils,graphics,grDevices,knitr
KNIT     := Rscript $(PKGS) -e "knit2html(\"$(SRC).Rmd\")"
PURL     := Rscript $(PKGS) -e "purl(\"$(SRC).Rmd\")"
VIEWHTML := open -a Safari
VIEWPDF  := open
EDIT     := aquamacs

all: $(SRC).html $(SRC).R openhtml

$(SRC).html: $(SRC).Rmd
	rm -f -r figure
	$(KNIT)
	rm -f .#confInt.Rmd

$(SRC).R: $(SRC).Rmd
	$(PURL)

openhtml:
	$(VIEWHTML) $(SRC).html
	$(EDIT) $(SRC).Rmd

SRC    := expDecay
KNIT   := Rscript -e "library(knitr); knit2pdf(\"$(SRC).Rnw\")"
PURL   := Rscript -e "library(knitr); purl(\"$(SRC).Rnw\")"
VIEW   := open
EDIT   := aquamacs

all: $(SRC).pdf $(SRC).R see edit clean

$(SRC).pdf: $(SRC).Rnw
	$(KNIT)

$(SRC).R: $(SRC).Rnw
	$(PURL)

see:
	$(VIEW) $(SRC).pdf

edit:
	$(EDIT) $(SRC).Rnw

clean:
	rm -r -f *.aux *.rel *.log *.toc

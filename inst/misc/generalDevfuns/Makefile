SRC    := generalDevfuns
KNIT   := Rscript -e "library(knitr); knit2pdf(\"$(SRC).Rnw\")"
PURL   := Rscript -e "library(knitr); purl(\"$(SRC).Rnw\")"
VIEW   := open
EDIT   := aquamacs

all: $(SRC).pdf $(SRC).R open clean

$(SRC).pdf: $(SRC).Rnw
	$(KNIT)

$(SRC).R: $(SRC).Rnw
	$(PURL)

open:
	$(VIEW) $(SRC).pdf
	$(EDIT) $(SRC).Rnw

clean:
	rm -r -f *.aux *.rel *.log *.toc


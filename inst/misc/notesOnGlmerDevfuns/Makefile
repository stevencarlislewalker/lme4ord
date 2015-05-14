SRC    := notesOnGlmerDevfuns
KNIT   := Rscript -e "library(\"knitr\"); knit2html(\"$(SRC).Rmd\")"
VIEW   := open -a Safari
EDIT   := aquamacs

all: $(SRC).html open

$(SRC).html: $(SRC).Rmd
	$(KNIT)

open:
	$(VIEW) $(SRC).html
	$(EDIT) $(SRC).Rmd

TEX = pdflatex
BIB = bibtex
SRC = modelDescription
BIBSRC = msLVM
VIEW = open -a Preview
EDIT = aquamacs
VER = 

all: $(SRC).pdf clean open

$(SRC).pdf:
	$(TEX) $(SRC)$(VER)
	$(TEX) $(SRC)$(VER)
	$(TEX) $(SRC)$(VER)

open:
	$(VIEW) $(SRC)$(VER).pdf
	$(EDIT) $(SRC)$(VER).tex

clean:
	rm -r -f *.aux *.log *.bbl *.blg

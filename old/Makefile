all: ph.pdf

figures := figure1.pdf plhalfar.pdf

figure1.pdf: plhalfar.py figure1.py
	./figure1.py

plhalfar.pdf: plhalfar.py
	./plhalfar.py

# TO BUILD ph.pdf:
# ice_bib.bib needs to be a link to the same object in pism-dev/doc/
# also need links from files in
#    http://www.igsoc.org/production/igs-v2_00-distrib.zip
# namely
#    igs.cls
#    igs.bst
#    igsnatbib.sty

ph.pdf: ph.aux ph.bbl ph.tex $(figures)
	pdflatex ph

ph.aux: ph.tex $(figures)
	pdflatex ph

ph.bbl: ph.aux ice_bib.bib
	bibtex ph

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *~ $(figures)

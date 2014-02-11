all: exmarine.pdf

# TO BUILD exmarine.pdf:
# ice_bib.bib needs to be a link to the same object in pism-dev/doc/
# also need links from files in
#    http://www.igsoc.org/production/igs-v2_00-distrib.zip
# namely
#    igs.cls
#    igs.bst
#    igsnatbib.sty

figures := twoparabolas.pdf bodverifbetaB.pdf bodverifthickvel.pdf verifN.pdf

twoparabolas.pdf: twoparabolas.py
	./twoparabolas.py

bodverifthickvel.pdf bodverifbetaB.pdf: exactfigures.py
	./exactfigures.py

verifN.pdf: conv.txt verifNfigure.py
	./verifNfigure.py conv.txt verifN.pdf

exmarine.pdf: exmarine.aux exmarine.bbl exmarine.tex $(figures)
	pdflatex exmarine

exmarine.aux: exmarine.tex $(figures)
	pdflatex exmarine

exmarine.bbl: exmarine.aux ice_bib.bib
	bibtex exmarine

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *~ $(figures)

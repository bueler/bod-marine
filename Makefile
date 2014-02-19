all: exmarine.pdf

# TO BUILD exmarine.pdf:
# ice_bib.bib needs to be a link to the same object in pism-dev/doc/
# also need links from files in
#    http://www.igsoc.org/production/igs-v2_00-distrib.zip
# namely
#    igs.cls
#    igs.bst
#    igsnatbib.sty

figures := twoparabolas.pdf verifN.pdf exactmarine-geometry.pdf

twoparabolas.pdf: twoparabolas.py
	./twoparabolas.py

exactmarine-geometry.pdf: marineshoot.py
	./marineshoot.py --saveroot exactmarine --noshoot
	# also generates exactmarine-M-B.pdf, exactmarine-beta.pdf

verifN.pdf: conv.txt verifNfigure.py
	./verifNfigure.py conv.txt verifN.pdf

exmarine.pdf: exmarine.aux exmarine.bbl exmarine.tex $(figures)
	pdflatex exmarine

exmarine.aux: exmarine.tex $(figures)
	pdflatex exmarine
	bibtex exmarine

exmarine.bbl: exmarine.aux ice-bib.bib
	bibtex exmarine

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~ $(figures) unnamed-*.pdf

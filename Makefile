all: exmarine.pdf

otherexmfigs = exactmarine-M-B.pdf exactmarine-beta-T.pdf \
          exactmarine-geometry-detail.pdf exactmarine-error.pdf \
          exactmarine-dt-adaptive.pdf

figures = twoparabolas.pdf convmarine.pdf exactmarine-geometry.pdf exactmarine-stiffness-ratio.pdf

twoparabolas.pdf: twoparabolas.py
	./twoparabolas.py

# also generates $(otherexmfigs)
exactmarine-geometry.pdf : marineshoot.py
	./marineshoot.py --saveroot exactmarine --figures

exactmarine-stiffness-ratio.pdf : linearization.py
	./linearization.py

convmarine.pdf: convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh
	./convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh

# TO BUILD exmarine.pdf:
# (1) ice_bib.bib needs to be a link to the same object in pism-dev/doc/
# (2) also need links from files in
#       http://www.igsoc.org/production/igs-v2_00-distrib.zip
#     namely  igs.cls igs.bst igsnatbib.sty
exmarine.pdf: exmarine.aux exmarine.bbl exmarine.tex $(figures)
	pdflatex exmarine

exmarine.aux: exmarine.tex $(figures)
	pdflatex exmarine
	bibtex exmarine

exmarine.bbl: exmarine.aux ice-bib.bib
	bibtex exmarine

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~ \
	$(figures) $(otherexmfigs) unnamed-*.pdf

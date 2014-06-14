all: exmarine.pdf

otheremfigs = em-M-B.pdf em-beta-T.pdf em-geometry-detail.pdf em-error.pdf \
          em-dt-adaptive.pdf

figures = twoparabolas.pdf convmarine.pdf em-geometry.pdf em-stiffness-ratio.pdf

twoparabolas.pdf: twoparabolas.py
	./twoparabolas.py

# also generates $(otheremfigs)
em-geometry.pdf : marineshoot.py
	./marineshoot.py --saveroot em --figures

em-stiffness-ratio.pdf : linearization.py
	./linearization.py

convmarine.pdf: convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh
	./convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh

# TO BUILD exmarine.pdf:
# (1) ice_bib.bib needs to be a link to the same object in pism-dev/doc/
# (2) also need links from files in
#       http://www.igsoc.org/production/igs-v2_00-distrib.zip
#     namely  igs.cls igs.bst igsnatbib.sty
#exmarine.pdf: exmarine.aux exmarine.bbl exmarine.tex $(figures)
exmarine.pdf: exmarine.aux exmarine.tex $(figures)
	pdflatex exmarine

exmarine.aux: exmarine.tex $(figures)
	pdflatex exmarine
#	bibtex exmarine

#exmarine.bbl: exmarine.aux ice-bib.bib
#	bibtex exmarine

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~ \
	$(figures) $(otheremfigs) unnamed-*.pdf

all: exmarine.pdf

otheremfigs = em-M-B.pdf em-beta-T.pdf em-geometry-detail.pdf em-error.pdf \
          em-dt-adaptive.pdf

figures = twoparabolas.pdf convmarine.pdf em-geometry.pdf em-stiffness-ratio.pdf

twoparabolas.pdf: twoparabolas.py
	./twoparabolas.py
	cp twoparabolas.pdf figure1.pdf

# also generates $(otheremfigs)
em-geometry.pdf : marineshoot.py
	./marineshoot.py --saveroot em --figures
	cp em-geometry.pdf figure2.pdf
	cp em-M-B.pdf figure3.pdf
	cp em-beta-T.pdf figure4.pdf
	cp em-geometry-detail.pdf figure5.pdf
	cp em-error.pdf figure6.pdf
	cp em-dt-adaptive.pdf figure7.pdf

em-stiffness-ratio.pdf : linearization.py
	./linearization.py
	cp em-stiffness-ratio.pdf figure8.pdf

convmarine.pdf: convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh
	./convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh
	cp convmarine.pdf figure9.pdf

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
	$(figures) $(otheremfigs) unnamed-*.pdf figure?.pdf

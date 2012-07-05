all: shelf_num.pdf

figures := testJ_H_contour.png testJ_nu_contour.png testJPISMaverr.png testJ_uv_quiver.png

# ice_bib.bib needs to be a link to the same object in pism-dev/doc/
shelf_num.pdf: shelf_num.aux shelf_num.bbl shelf_num.tex $(figures)
	pdflatex shelf_num

shelf_num.aux: shelf_num.tex $(figures)
	pdflatex shelf_num

shelf_num.bbl: shelf_num.aux ice_bib.bib
	bibtex shelf_num

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg


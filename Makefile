all: ph.pdf


figures := figure1.pdf plhalfar.pdf twoparabolas.pdf bodverifbetaB.pdf \
           bodverifthickvel.pdf verifN.pdf

figure1.pdf: plhalfar.py figure1.py
	./figure1.py

plhalfar.pdf: plhalfar.py
	./plhalfar.py

twoparabolas.pdf: plhalfar.py twoparabolas.py
	./twoparabolas.py

bodverifthickvel.pdf bodverifbetaB.pdf: bodverifpicture.py
	./bodverifpicture.py

verifN.pdf: conv.txt verifNfigure.py
	./verifNfigure.py conv.txt verifN.pdf

# ice_bib.bib needs to be a link to the same object in pism-dev/doc/
ph.pdf: ph.aux ph.bbl ph.tex $(figures)
	pdflatex ph

ph.aux: ph.tex $(figures)
	pdflatex ph

ph.bbl: ph.aux ice_bib.bib
	bibtex ph

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *~ $(figures)

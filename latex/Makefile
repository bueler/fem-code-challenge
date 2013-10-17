# Generic make file for LaTeX: requires GNU make

all: obstacleDOC.pdf

obstacleDOC.pdf: obstacleDOC.tex icefe.bib
	pdflatex $<
	bibtex $(<:.tex=.aux)
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~
	@rm -f *converted-to.pdf


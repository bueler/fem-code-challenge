# Generic make file for LaTeX: requires GNU make

all: $(TEXFILE:.tex=.pdf)

%.pdf: %.tex icefe.bib
	pdflatex $<
	bibtex $(<:.tex=.aux)
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~


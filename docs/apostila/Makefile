# Apostila do GHOST. Requer XeLaTeX (as fontes da identidade estao em fontes/).
APOSTILA = apostila

all: $(APOSTILA).pdf

$(APOSTILA).pdf: $(APOSTILA).tex preambulo.tex capa.tex $(wildcard cap*.tex)
	xelatex -interaction=nonstopmode $(APOSTILA).tex
	xelatex -interaction=nonstopmode $(APOSTILA).tex

ver: $(APOSTILA).pdf
	xdg-open $(APOSTILA).pdf

clean:
	rm -f *.aux *.log *.out *.toc *.fls *.fdb_latexmk *.synctex.gz
	rm -f *.png teste.* a.log t.log

.PHONY: all ver clean

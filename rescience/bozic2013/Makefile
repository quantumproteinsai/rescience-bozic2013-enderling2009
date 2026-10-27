# Build the article with XeLaTeX (ReScience C template). metadata.tex is written by hand
# from metadata.yaml (same definitions as the template's yaml-to-latex.py produces).
article.pdf: article.tex content.tex bibliography.bib metadata.tex
	latexmk -pdf -pdflatex="xelatex -interaction=nonstopmode" article.tex

clean:
	latexmk -CA

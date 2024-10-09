# Makefile for LaTeX

# Variables
TEX = pdflatex
FILE = zuoye1

all: $(FILE).pdf

# Compile .tex to .pdf
$(FILE).pdf: $(FILE).tex
	$(TEX) $(FILE).tex

# Clean up auxiliary files
clean:
	rm -f *.aux *.log *.out *.toc

# Clean up everything including the PDF
cleanall: clean
	rm -f $(FILE).pdf

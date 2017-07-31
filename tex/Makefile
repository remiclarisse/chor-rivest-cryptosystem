# Variables
LOGBOOK=brouillon
RAPPORT=rapport-stage

# Special rules and targets
.PHONY: all clean help open

# Rules and targets
all: $(RAPPORT)

$(LOGBOOK): clean
	pdflatex $(LOGBOOK).tex
	bibtex $(LOGBOOK).aux
	pdflatex $(LOGBOOK).tex
	pdflatex $(LOGBOOK).tex

$(RAPPORT): clean
	pdflatex $(RAPPORT).tex
	bibtex $(RAPPORT).aux
	pdflatex $(RAPPORT).tex
	pdflatex $(RAPPORT).tex

open: $(RAPPORT)
	@evince $(RAPPORT).pdf&

clean:
	@echo "./tex/ : Cleaning..."
	@rm -f *~ *.aux *.bbl *.blg *.toc *.log *.gz *.out *.pdf *.gz\(busy\)

help:
	@echo "all:\tRun the whole build to make the PDF file '$(RAPPORT).pdf'."
	@echo "open:\tBuild the PDF file '$(RAPPORT).pdf' and open it."
	@echo "clean:\tRemove all files produced by the compilation."
	@echo "help:\tDisplay the main targets of this Makefile with a short despriction."

# Variables
LOGBOOK=chor-rivest
RAPPORT=rapport-stage

# Special rules and targets
.PHONY: all clean help open

# Rules and targets
all: clean $(LOGBOOK) $(RAPPORT)

$(LOGBOOK):
	pdflatex $(LOGBOOK).tex > /dev/null
	bibtex $(LOGBOOK).aux > /dev/null
	pdflatex $(LOGBOOK).tex > /dev/null
	pdflatex $(LOGBOOK).tex > /dev/null

$(RAPPORT):
	pdflatex $(RAPPORT).tex > /dev/null
	bibtex $(RAPPORT).aux > /dev/null
	pdflatex $(RAPPORT).tex > /dev/null
	pdflatex $(RAPPORT).tex > /dev/null

open: $(LOGBOOK) $(RAPPORT)
	@evince $(LOGBOOK).pdf&
	@evince $(RAPPORT).pdf&

clean:
	@echo "./tex/ : Cleaning..."
	@rm -f *~ $(LOGBOOK).aux $(LOGBOOK).bbl $(LOGBOOK).blg $(LOGBOOK).toc *.log *.gz $(LOGBOOK).out
	@rm -f *~ $(RAPPORT).aux $(RAPPORT).bbl $(RAPPORT).blg $(RAPPORT).toc *.log *.gz $(RAPPORT).out

help:
	@echo "all:\tRun the whole build to make the PDF files."
	@echo "open:\tBuild the PDF files and open it."
	@echo "clean:\tRemove all files produced by the compilation."
	@echo "help:\tDisplay the main targets of this Makefile with a short despriction."

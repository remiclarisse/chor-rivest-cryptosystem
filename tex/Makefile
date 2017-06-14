# Variables
FILE=chor-rivest

# Special rules and targets
.PHONY: all clean help open

# Rules and targets
all: $(FILE)

$(FILE): clean
	pdflatex $(FILE).tex > /dev/null
	bibtex $(FILE).aux > /dev/null
	pdflatex $(FILE).tex > /dev/null
	pdflatex $(FILE).tex > /dev/null

open: $(FILE)
	@evince $(FILE).pdf&

clean:
	@echo "./tex/ : Cleaning..."
	@rm -f *~ $(FILE).aux $(FILE).bbl $(FILE).blg $(FILE).toc *.log *.gz

help:
	@echo "all:\tRun the whole build to make the PDF file."
	@echo "open:\tBuild the PDF file and open it."
	@echo "clean:\tRemove all files produced by the compilation."
	@echo "help:\tDisplay the main targets of this Makefile with a short despriction."

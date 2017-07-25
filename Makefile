# Variables
EXE=rho-pollard
FILE=rapport-stage
DLP=computing-dlp
IMPLEMCHORRIVEST=implementation-chor-rivest
VAUDENAY=with-vaudenay-attack
JOUX=with-joux-algorithm

# Special rules and targets
.PHONY: all clean help info open

# Rules and tagets
all: clean
	@mkdir $(FILE)
	@cd tex && $(MAKE) $(FILE)
	@cp -f tex/$(FILE).pdf ./$(FILE)/
	@mkdir $(FILE)/$(IMPLEMCHORRIVEST)
	@mkdir $(FILE)/$(IMPLEMCHORRIVEST)/$(VAUDENAY)
	@cp -f sage/chor-rivest-vaudenay/* ./$(FILE)/$(IMPLEMCHORRIVEST)/$(VAUDENAY)
	@mkdir $(FILE)/$(IMPLEMCHORRIVEST)/$(JOUX)
	@cp -f sage/chor-rivest-joux/* ./$(FILE)/$(IMPLEMCHORRIVEST)/$(JOUX)
	@mkdir $(FILE)/$(DLP)
	@cp -f sage/DLP/hellman-reyneri.sage ./$(FILE)/$(DLP)/
	@cp -f sage/DLP/pohlig-hellman.sage ./$(FILE)/$(DLP)/
	@cp -f sage/DLP/README.md ./$(FILE)/$(DLP)/
	@cd c && $(MAKE) nodebug
	@mkdir $(FILE)/$(EXE)
	@cp -f c/$(EXE) ./$(FILE)/$(EXE)/
	@cp -f c/README.md ./$(FILE)/$(EXE)/
	@cp -f c/test/primes ./$(FILE)/$(EXE)/

$(EXE):
	@cd c && $(MAKE) nodebug
	@cp -f c/$(EXE) ./

pdf:
	@cd tex && $(MAKE) $(FILE)
	@cp -f tex/$(FILE).pdf ./

clean:
	@cd c && $(MAKE) clean
	@cd tex && $(MAKE) clean
	@echo "./ : Cleaning..."
	@rm -f *~ $(EXE) *.pdf *.log
	@rm -rf ./$(FILE)

open: pdf
	@evince $(FILE).pdf&

info:
	@more README.md

help:
	@echo "all:\t\tRun the target build."
	@echo "$(EXE):\tBuid executable file $(EXE) recursively in the directory."
	@echo "pdf:\t\tBuid PDF text file $(FILE).pdf recursively in the directory."
	@echo "build:\t\tBuild all executable file(s) and PDF text file(s) recursively in the directory."
	@echo "clean:\t\tRemove all files produced by the compilation."
	@echo "info:\t\tGive info about this project."
	@echo "help:\t\tDisplay the main targets of this Makefile with a short despriction."

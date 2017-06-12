# Variables
EXE=rho-pollard
FILE=chor-rivest

# Special rules and targets
.PHONY: all clean help info

# Rules and tagets
all: $(EXE) $(FILE)

$(EXE):
	@cd c && $(MAKE)
	@cp -f c/$(EXE) ./

$(FILE):
	@cd tex && $(MAKE)
	@cp -f tex/$(FILE).pdf ./

clean:
	@cd c && $(MAKE) clean
	@cd tex && $(MAKE) clean
	@echo "./ : Cleaning..."
	@rm -f *~ $(EXE) $(FILE).pdf

info:
	@more README.md

help:
	@echo "all:\t\tRun the target build."
	@echo "$(EXE):\tBuid executable file $(EXE) recursively in the directory."
	@echo "$(FILE):\tBuid PDF text file $(FILE).pdf recursively in the directory."
	@echo "build:\t\tBuild all executable file(s) and PDF text file(s) recursively in the directory."
	@echo "clean:\t\tRemove all files produced by the compilation."
	@echo "info:\t\tGive info about this project."
	@echo "help:\t\tDisplay the main targets of this Makefile with a short despriction."

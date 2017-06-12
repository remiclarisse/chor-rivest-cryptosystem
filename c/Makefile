# Variables
EXE=rho-pollard

# Special rules and targets
.PHONY: all build clean help nodebug

# Rules and tagets
all: build

build:
	@cd src && $(MAKE)
	@cp -f src/$(EXE) ./

clean:
	@cd src && $(MAKE) clean
	@echo "./c/ : Cleaning..."
	@rm -f $(EXE) $(EXE)_nodebug

nodebug:
	@cd src && $(MAKE) nodebug
	@cp -f src/$(EXE)_nodebug ./

help:
	@echo "all:\tRun the target build."
	@echo "build:\tBuild the executable file rho-pollard recursively in the directory"
	@echo "clean:\tRemove all files produced by the compilation."
	@echo "help:\tDisplay the main targets of this Makefile with a short despriction."

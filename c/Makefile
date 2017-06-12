# Variables
EXE=rho-pollard

# Special rules and targets
.PHONY: all build clean help

# Rules and tagets
all: build

build:
	@cd src && $(MAKE)
	@cp -f src/$(EXE) ./

clean:
	@cd src && $(MAKE) clean
	@echo "./c/ : Cleaning..."
	@rm -f $(EXE)

help:
	@echo "all:\t Run the target build."
	@echo "build:\t Build the executable file rho-pollard recursively in the directory"
	@echo "clean:\t Remove all files produced by the compilation."
	@echo "help:\t Display the main targets of this Makefile with a short despriction."

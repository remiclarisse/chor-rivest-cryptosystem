# Variables
EXE=rho-pollard

# Usual compilation flags
CFLAGS=-std=c11 -Wall -Wextra -g3
CPPFLAGS=-DDEBUG
LDFLAGS=-lm

# Special rules and targets
.PHONY: all clean help nodebug

# Rules and targets
all: $(EXE)

$(EXE): rho-pollard.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

rho-pollard.o: rho-pollard.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

nodebug:
	gcc -std=c11 -o $(EXE)_nodebug $(EXE).c

clean:
	@echo "./c/src/ : Cleaning..."
	@rm -f *~ *.o $(EXE) $(EXE)_nodebug

help:
	@echo "all:\t Run the whole build of pollard's rho logarithm computing program."
	@echo "reversi: Build the executable file from rho-pollard.c"
	@echo "clean:\t Remove all files produced by the compilation."
	@echo "help:\t Display the main targets of this Makefile with a short despriction."

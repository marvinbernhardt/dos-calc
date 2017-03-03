PREFIX = $$HOME
OBJ = dos-calc.c
PROGRAM = dos-calc

CC = /usr/bin/gcc
CFLAGS = -Wall -O3 -fopenmp -std=gnu99
LDFLAGS = -lgromacs -llapacke -lcblas -lblas -lfftw3f -lm
INC = -I $${GROMACS_DIR}/include
LIB = -L $${GMXLDLIB}

$(PROGRAM): $(OBJ) velocity-decomposition.c fft.c linear-algebra.c verbPrintf.c
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJ) $(INC) $(LIB) $(LDFLAGS)

.PHONY: install
install: $(PROGRAM)
	mkdir -p $(PREFIX)/bin
	cp $(PROGRAM) $(PREFIX)/bin/$(PROGRAM)

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/$(PROGRAM)

.PHONY: clean
clean:
	rm -f $(PROGRAM)

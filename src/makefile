PREFIX := $$HOME
OBJ := dos-calc.c
PROGRAM := dos-calc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

CFLAGS := -Wall -O3 -fopenmp -std=gnu99 -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS := -lgromacs -llapacke -llapack -lcblas -lblas -lfftw3f -lm -lgfortran
INC := -I $${GROMACS_DIR}/include
LIB := -L $${GMXLDLIB}

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

#/bin/bash

gcc dos-calc-wrapper.c \
    -Wall \
    -march=nehalem \
    -O3 \
    -fopenmp \
    -std=gnu99 \
    -o $HOME/bin/dos-calc-openmp \
    -I $HOME/offline/software/gromacs-2016.1/include \
    -L $GMXLDLIB \
    -lgromacs -llapacke -lcblas -lblas -lfftw3f -lm 

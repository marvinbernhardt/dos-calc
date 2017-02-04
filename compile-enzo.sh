#!/bin/bash

gcc dos-calc-wrapper.c \
    -Wall \
    -std=gnu99 \
    -o $HOME/bin/dos-calc \
    -I $HOME/software/gromacs-2016.1/include \
    -L $GMXLDLIB \
    -L $HOME/software/lib \
    -I $HOME/software/include \
    -lgromacs -llapacke -llapack -lcblas -lblas -lfftw3f -lm -lgfortran

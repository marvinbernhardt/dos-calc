#/bin/bash

gcc dos-calc-wrapper.c \
    -Wall \
    -O3 \
    -std=gnu99 \
    -o $HOME/bin/dos-calc-devel \
    -I $HOME/offline/software/gromacs-2016.1/include \
    -L $GMXLDLIB \
    -lgromacs -llapacke -lcblas -lblas -lfftw3f -lm 

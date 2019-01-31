#!/bin/bash
# This script will install dos-calc in $HOME/bin/dos-calc
set -euo pipefail

mkdir -p build
pushd build
    rm -rf *
    cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME ..
    make
    make install
popd

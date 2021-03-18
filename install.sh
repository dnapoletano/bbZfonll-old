#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/gfortran/lib:/usr/local/lib:/usr/lib64
rm -r lib
mkdir lib
cd src
make clean
make
make python
cd -

cp src/_bbHFONLL.so src/bbHFONLL.py lib/.

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:lib

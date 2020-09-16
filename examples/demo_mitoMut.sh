#!/bin/sh

cd ../

## Input files
vireoMT -m data/mitoDNA/cellSNP.tag.AD.mtx,data/mitoDNA/cellSNP.tag.DP.mtx \
    -o data/mitoMutOUT -p 4
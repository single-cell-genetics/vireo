#!/bin/sh

## MODE 1: no genotype
CELL_FILE=data/cells.cellSNP.vcf.gz
OUT_DIR=data/cellSNP_noGT
vireo -c $CELL_FILE -N 4 -o $OUT_DIR --randSeed 2

# CELL_FILE=data/cells.donorid.vcf.gz
# OUT_DIR=data/donorid_noGT
# vireo -c $CELL_FILE -N 3 -o $OUT_DIR

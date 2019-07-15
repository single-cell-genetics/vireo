#!/bin/sh

## Input files
CELL_FILE=data/cells.cellSNP.vcf.gz
DONOR_FILE=data/donors.cellSNP.vcf.gz

mkdir data/outs/

## MODE 1: no genotype
OUT_DIR=data/outs/cellSNP_noGT
vireo -c $CELL_FILE -N 4 -o $OUT_DIR --randSeed 2 #--amplifyK 1.5

## MODE 2: given genotype
OUT_DIR=data/outs/cellSNP_PL
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR --randSeed 2 #--genoTag PL

## MODE 3: given genotype with learn
OUT_DIR=data/outs/cellSNP_learn
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR --randSeed 2 -N 4 --forceLearnGT


# CELL_DIR=data/outs/sparseMat
# OUT_DIR=data/outs/cellSNP_PL2
# vireo -c $CELL_DIR -d $DONOR_FILE -o $OUT_DIR --randSeed 2


# CELL_FILE=data/cells.donorid.vcf.gz
# OUT_DIR=data/donorid_noGT
# vireo -c $CELL_FILE -N 3 -o $OUT_DIR

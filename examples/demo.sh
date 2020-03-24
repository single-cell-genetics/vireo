#!/bin/sh

cd ../

## Input files
CELL_DIR=data/cellSNP_mat
CELL_FILE=data/cells.cellSNP.vcf.gz
DONOR_FILE=data/donors.cellSNP.vcf.gz
DONOR_FILE_PART=data/donors.two.cellSNP.vcf.gz

mkdir data/outs/

## MODE 1: no donor genotype
OUT_DIR=data/outs/cellSNP_noGT
vireo -c $CELL_DIR -N 4 -o $OUT_DIR --randSeed 2 #--ASEmode # --extraDonor 0
# vireo -c $CELL_DIR -N 4 -o $OUT_DIR --randSeed 2 --extraDonor 1 #--ASEmode

## MODE 2: given donor genotype
OUT_DIR=data/outs/cellSNP_PL
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR -N 4 --randSeed 2 #--genoTag PL


## MODE 3: given part donor genotype
OUT_DIR=data/outs/cellSNP_part
vireo -c $CELL_FILE -d $DONOR_FILE_PART -o $OUT_DIR -N 4 --randSeed 2


## MODE 4: given donor genotype but not perfect, use as prior to learn
OUT_DIR=data/outs/cellSNP_learn
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR --randSeed 2 -N 4 --forceLearnGT


## MODE 5: given donor genotype but too many
OUT_DIR=data/outs/cellSNP_PL3
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR -N 3 --randSeed 2 #--genoTag PL


## Generating genotype barcodes
donor_vcf=data/outs/cellSNP_noGT/GT_donors.vireo.vcf.gz
GTbarcode -i $donor_vcf -o data/outs/cellSNP_noGT/GT_barcodes.tsv --randSeed 1 # --noHomoAlt

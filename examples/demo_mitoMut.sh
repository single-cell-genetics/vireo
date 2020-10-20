#!/bin/sh

#cd ../

## Input files
TEST_SPARSE_AD=data/mitoDNA/cellSNP.tag.AD.mtx
TEST_SPARSE_DP=data/mitoDNA/cellSNP.tag.DP.mtx
TEST_SPARSE_VCF=data/mitoDNA/kim_cellSNP.cells.vcf.gz

#test input AD and DP sparse matrices 
vireoMT -m $TEST_SPARSE_AD,$TEST_SPARSE_DP -o data/mitoDNA/mitoMutOUT_sparse -p 5

vireoMT --vcfData $TEST_SPARSE_VCF -o data/mitoDNA/mitoMutOUT_vcf -p 8
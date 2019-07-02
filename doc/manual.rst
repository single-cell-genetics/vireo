======
Manual
======

Demultiplexing requires two read or UMI count matrices for each variant in each
cell: `A` for alternative allele and `D` depth (i.e., summary of alternative 
and reference alleles).


In many biomedical studies, biological replications are important to eliminate 
unwanted variations, e.g., genetic variations, and often the genetic variation 
is not the main focus, and genotype reference is probably not available.

For Vireo, known genotype reference is not necessary to enjoy the benefits of 
multiplexed scRNA-seq. The only requirement is a cell VCF file and the number 
of pooled donors (actually Vireo can detect the latter automatically).

  * cell VCF file: a variant call format (VCF) file from which we can extract 
    the variant x cell matrices of integer count of the number of reads 
    supporting the alternative allele for each variant in each cell and the 
    total number of reads overlapping each variant in each cell. This file can 
    be piled up from bam file on a list of common SNPs by the multifaceted 
    bcftools or our tailored designed cellSNP, Python package based on pysam.

Here, we demonstrate the use of Vireo to assign 384 cells to 3 donors by 2,171 
SNPs, which can be loaded from cardelino package directly, via a function 
load_cellSNP_vcf based on vcfR package.


referece allele in each

Once the genotypes for each cell have been obtained (see genotype_ page), 
specifically the read or UMI counts for each variant (both alternative and 
referece allele in each we 
can demultiplex

Demultiplexing from allelic expression
--------------------------------------

This python package offers a set of utilities functions and an executable 
command line `vireo` for donor deconvolution in any of these four situations:

1) without any genotype: 

   ::

      vireo -c $CELL_FILE -N $n_donor -o $OUT_DIR

2) with genotype for all samples (GT, GP, or PL)

   ::

      vireo -c $CELL_FILE -d $DONOR_GT_FILE -o $OUT_DIR

3) with genotype for part of the samples

   ::

      vireo -c $CELL_FILE -d $DONOR_GT_FILE -o $OUT_DIR -N $n_donor 

4) with genotype but not confident

   ::

      vireo -c $CELL_FILE -d $DONOR_GT_FILE -o $OUT_DIR --forceLearnGT

For details, see the full manual_ or type "vireo -h" for all arguments. We also 
provide a demo.sh_ for running the test data sets in this repo.

.. _manual: https://vireoSNP.readthedocs.io/en/latest/manual.html
.. _demo.sh: https://github.com/huangyh09/vireo/blob/master/demo.sh
.. _genotype: https://vireoSNP.readthedocs.io/en/latest/genotype.html


======
Manual
======

Demultiplexing requires two count matrices (variant-by-cell) of reads or UMIs 
for each variant in each cell: ``A`` for alternative allele and ``D`` depth 
(i.e., summary of alternative and reference alleles). These two matrices can be 
obtained by genotyping a list of variants in each cell. We provide a guideline 
for cellular genotyping_ with a recommendation of cellSNP_ that is developed by 
us, too.

Once the genotypes for each cell have been obtained, e.g., in VCF format, or two
sparse matrices ``A`` and ``D``, we can apply Vireo for demultiplexing.


Demultiplexing single cells
===========================

By default, Vireo works without any known genotype information for pooled 
samples. However, if any of genotype of these samples are known or can be 
obtained, e.g., by bulk RNA-seq, exome-seq, it is still useful to add them, not
only allowing us to align the deconvoluted samples to its identity, but also can 
benefits the doublets identification, especially if the coverage or the loaded 
cells per sample is low.

Depending the availability of genotype information, we provide four strategies 
to demultiplex scRNA-seq data.

1) without any genotype: 

   ::

      vireo -c $CELL_FILE -N $n_donor -o $OUT_DIR

2) with genotype for all samples (GT, GP, or PL)

   ::

      vireo -c $CELL_FILE -d $DONOR_GT_FILE -o $OUT_DIR

3) with genotype for part of the samples

   ::

      vireo -c $CELL_FILE -d $DONOR_GT_FILE -o $OUT_DIR -N $n_donor 

4) with genotype but not confident (or only for subset of SNPs)

   ::

      vireo -c $CELL_FILE -d $DONOR_GT_FILE -o $OUT_DIR --forceLearnGT


All arguments
=============

Type ``vireo -h`` for details of all arguments:

.. code-block:: html

   Usage: vireo [options]

   Options:
      -h, --help            show this help message and exit
      -c CELL_FILE, --cellFile=CELL_FILE
                              The cell genotype file in VCF format
      -N N_DONOR, --nDonor=N_DONOR
                              Number of donors to demultiplex; can be larger than
                              provided in donor_file
      -o OUT_DIR, --outDir=OUT_DIR
                              Dirtectory for output files [default:
                              $cellFilePath/vireo]

      Optional arguments:
      --noDoublet         If use, not checking doublets.
      -d DONOR_FILE, --donorFile=DONOR_FILE
                              The donor genotype file in VCF format. Please filter
                              the sample and region with bcftools -s and -R first!
      -t GENO_TAG, --genoTag=GENO_TAG
                              The tag for donor genotype: GT, GP, PL [default: PL]
      -M N_INIT, --nInit=N_INIT
                              Number of random initializations [default: 2 (GT) or
                              50 (no GT)]
      --amplifyK=K_AMPLIFY
                              Pre-cluster with amplified K [default: 1.0 (GT) or 1.2
                              (no GT)]
      --forceLearnGT      If use, treat donor GT as prior only.
      --randSeed=RAND_SEED
                              Seed for random initialization [default: none]


Example data
============

In order to test vireo and illustrate the usage, we provide a test `data set`_,
also some `demo scripts`_.

This example data set contains 952 cells from 4 samples. The genotypes for these
four samples are also provided.

.. _genotyping: https://vireoSNP.readthedocs.io/en/latest/genotype.html
.. _cellSNP: https://github.com/huangyh09/cellSNP
.. _demo scripts: https://github.com/huangyh09/vireo/blob/master/demo.sh
.. _data set: https://github.com/huangyh09/vireo/tree/master/data

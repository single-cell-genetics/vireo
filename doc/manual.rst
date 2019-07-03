======
Manual
======

Demultiplexing requires two count matrices of reads or UMIs for each variant in 
each cell: `A` for alternative allele and `D` depth (i.e., summary of alternative 
and reference alleles). These two matrices can be obtained by genotyping a list 
of variants in each cell. We provide a guideline for genotyping here_ with a 
recommendation to cellSNP_ that is also developed by us.

Once the genotypes for each cell have been obtained, e.g., in VCF format, or two
sparse matrices `A` and `D`, we can apply Vireo for demultiplexing.

By default, Vireo works without any known genotype information for pooled 
samples. However, if any of genotype of these samples are known or can be 
obtained, e.g., by bulk RNA-seq, exome-seq, it is still useful to add them, not
only allowing us to align the deconvoluted samples to its identity, but also can 
benefits the doublets identification, especially if the coverage or the loaded 
cells per sample is low.

Depending the availability of genotype information, we provide four strategies 
to demultiplex scRNA-seq data.

Demultiplexing from allelic expression
--------------------------------------

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

For details, type "vireo -h" for all arguments. We also provide a demo.sh_ for 
running the test `data sets`_.

.. _here: https://vireoSNP.readthedocs.io/en/latest/genotype.html
.. _cellSNP: https://github.com/huangyh09/cellSNP
.. _demo.sh: https://github.com/huangyh09/vireo/blob/master/demo.sh
.. _data sets: https://github.com/huangyh09/vireo/tree/master/data

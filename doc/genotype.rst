==========
Genotyping
==========

Genotyping (or piling up) a list of common variants on each cell is a pre-step 
for demultiplexing them with Vireo. This step requires some bioinformatics 
efforts, but thanks to many developers in this community, there are a few 
good existing software to use.

Genotyping cells can be divided into the following two sub-steps, and in 
different situations, the strategy may need to be customised for ensuring 
high quality genotyping data to demultiplex cells.

**Recommended strategies for genotyping cells:**

* For human or genotyped species: `variant list`_ (given) + cellsnp-lite_ (typing).
* For species without known common variants: cellsnp-lite_ (calling with mode 2b
  & typing with mode 1a). freebayes_ is an alternative to cellsnp-lite_ mode 2b
  for calling of heterozygous SNPs.

.. note::
   cellSNP_ was initially developed in Python based on pysam, which is 
   convenient for a few thousand cells but becomes the computational bottleneck 
   for large number of cells. Therefore, we have re-implemented it to 
   cellsnp-lite_ in C/C++ with ~5x faster and ~50x less memory.


1. Identify candidate SNPs
===========================
There are multiple ways to identify candidate SNPs, each has unique properties 
and may suit situations or species differently. Here, we listed two common 
strategies.

**Option 1): using a common variants**

The best way is to genotype the each of the pooled samples either by genotyping 
array or exome/genome/RNA sequencing (or other omics data), and use this list of 
distinguish variants to as candidate to genotype each cell. However, this can 
be costly, and not necessary in most cases.

For human, a very comprehensive list of common variants have been identified 
by international efforts, e.g., `10000 genome project`_ and gnomAD_ which gives 
a few millions common SNPs to genotype on each cell. The benefits include the 
reduced confounders, e.g., caused by RNA editing. We normally recommend this if 
for human, and we provide some `pre-processed SNP list`_.

.. note::
  Here are some tips if you have genotypes for input donors:
  
  1. Imputation can be helpfule if you obtained genotypes for each human  
     individual from SNP-array or Whole exome-seq. 
  2. Selection of informative SNPs is useful for filtering out SNPs with 
     identical or very similar genotypes in all mixed donors. You may consider
     `AC`, `AF` or similar tag in you donor VCF file. bcftools_ is a very 
     useful tool for such preprocessing.


**Option 2): Calling variants from scRNA-seq**

Other than human, most species may not have a well-defined common variants, 
hence the best way is to call the variants from pooled scRNA-seq directly.

The package freebayes_ is an often choice, which designed to find small 
polymorphisms, specifically SNPs, indels, MNPs, and complex events smaller than 
the length of a short-read sequencing alignment. Importantly, freebayes_ has 
a set of options to filter reads and variants.

We also recommend an alternative method cellsnp-lite_ that is developed by us.  
Its mode 2b has a similar feature to pileup the whole genome and identify the 
heterozygous variants in the pooled samples. It has highly comparable 
accuracry to freebayes_ and bcftools mpileup_ and achieves 5-10x speedups.


2. Genotype each cell
=====================

Once a list of candidate variants are found, it is more straightforward to 
genotype each cell. We provide three common methods, with recommendation to 
our cellsnp-lite_ to seamlessly with Vireo.

* The famous mpileup_ from bcftools / samtools is often a good choice. However, 
  this can be slow, as it doesn't allow parallel  computing and it doesn't 
  support the cell barcodes and UMI tag in the pooled BAM file for many cells.

* This limitation motivates us to develop cellSNP_, a pysam wrap (now 
  cellsnp-lite_ in C/C++) to pile up the 
  variants in each cell. The benefits include parallel computing, taking cell 
  barcoding tag and support UMIs.

* Alternatively, vartrix_ is also an option to genotype cells in 10x Genomics 
  data. 

Once the genotype (mainly pileup) has been achieved, it can be used to 
demultiplex the pooled cells, see the manual_.


.. _gnomAD: https://gnomad.broadinstitute.org/
.. _10000 genome project: http://www.internationalgenome.org/
.. _variant list: https://sourceforge.net/projects/cellsnp/files/SNPlist/
.. _pre-processed SNP list: https://sourceforge.net/projects/cellsnp/files/SNPlist/
.. _freebayes: https://github.com/ekg/freebayes
.. _cellSNP: https://github.com/single-cell-genetics/cellSNP
.. _cellsnp-lite: https://cellsnp-lite.readthedocs.io/en/latest/manual.html
.. _mpileup: http://www.htslib.org/doc/bcftools.html
.. _vartrix: https://github.com/10XGenomics/vartrix
.. _manual: https://vireosnp.readthedocs.io/en/latest/manual.html
.. _bcftools: http://samtools.github.io/bcftools/bcftools.html
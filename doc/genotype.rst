==========
Genotyping
==========

Genotyping (or piling up) a list of common variants on all cells is the first 
step for demultiplexing with Vireo. This step requires some bioinformatics 
efforts, but thanks to many developers in this community, there are a few 
good existing software to use.

Genotyping cells can be divided into the following two sub-steps, and in 
different situations, the strategy may need to be optimised for a 
high quality genotyping data to demultiplex cells.

**Recommended strategies for genotyping cells:**

* For human or genotyped species: `variant list`_ (given) + cellSNP_ (typing).
* For species without known common variants: freebayes_ (calling) + cellSNP_ 
  (typing)

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


**Option 2): Calling variants from scRNA-seq**

Other than human, most species may not have a well-defined common variants, 
hence the best way is to call the variants from pooled scRNA-seq directly.

The package freebayes_ is an often choice, which designed to find small 
polymorphisms, specifically SNPs, indels, MNPs, and complex events smaller than 
the length of a short-read sequencing alignment. Importantly, freebayes_ has 
a set of options to filter reads and variants.

Alternatively, cellSNP_ has a similar feature (still under development) to 
pileup the whole genome and identify the heterozygous variants in the pooled 
samples. However, this mode doesn't have a sophisticated model and filtering 
strategy for indentifying candidate SNPs, and may struggle with confounders 
from RNA editing.


2. Genotype each cell
=====================

Once a list of candidate variants are found, it is more straightforward to 
genotype each cell. We provide three common methods, with recommendation to 
``cellSNP``, which is developed by us.

* The famous mpileup_ from bcftools / samtools is often a good choice. However, 
  this can be slow, as it doesn't allow parallel  computing and it doesn't 
  support the cell barcodes and UMI tag in the pooled BAM file for many cells.

* This limitation motivates us to develop cellSNP_, a pysam wrap to pile up the 
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
.. _cellSNP: https://github.com/huangyh09/cellSNP
.. _mpileup: http://www.htslib.org/doc/bcftools.html
.. _vartrix: https://github.com/10XGenomics/vartrix
.. _manual: https://vireosnp.readthedocs.io/en/latest/manual.html
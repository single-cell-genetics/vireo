======================================================
vireo: donor deconvolution for pooled single-cell data
======================================================

Vireo: Variational Inference for Reconstructing Ensemble Origin by expressed 
SNPs in multiplexed scRNA-seq data. 

The name vireo follows the theme from cardelino_ (for clone deconvolution), 
while the Python package name is vireoSNP_ to aviod name confilict on PyPI.

.. _cardelino: https://github.com/PMBio/cardelino
.. _vireoSNP: https://pypi.org/project/vireoSNP


Installation
============

The easiest way is to install via PyPI_ by typing this line is terminal:

  ::

    pip install vireoSNP

Alternatively, you can always download this repository and install manually:

  ::

    python setup.py install

For more options of installation, see the full installation_.

.. _PyPI: https://pypi.org/project/vireoSNP
.. _manual: https://vireoSNP.readthedocs.io/en/latest/manual.html
.. _installation: https://vireoSNP.readthedocs.io/en/latest/install.html


Usage and manual
================

Genotyping for each cell (pre-step)
-----------------------------------
There might be some bioinformatics efforts in this step, however, a few existing 
software can provide a solution. There are often two steps for this:

1) identify candidate SNPs: `known common SNPs`_ / freebayes_ / cellSNP_
2) genotype candidate SNPs in each cell: cellSNP_ / vartrix_ / `bcftools mpileup`_

See more introduction in the genotyping_ section.

.. _known common SNPs: https://github.com/huangyh09/cellSNP#list-of-candidate-snps
.. _freebayes: https://github.com/ekg/freebayes
.. _cellSNP: https://github.com/huangyh09/cellSNP
.. _vartrix: https://github.com/10XGenomics/vartrix
.. _bcftools mpileup: http://www.htslib.org/doc/bcftools.html
.. _genotyping: https://vireoSNP.readthedocs.io/en/latest/genotype.html


Demultiplexing from allelic expression
--------------------------------------

This python package offers a set of utilities functions and an executable 
command line `vireo` for donor deconvolution in any of these four situations:

1) without any genotype: 

   ::

      vireo -c $CELL_FILE -N $n_donor -o $OUT_DIR

2) with genotype for all samples (tag via -t: GT, GP, or PL)

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


Reference
=========

Yuanhua Huang, Davis J. McCarthy, and Oliver Stegle. `Vireo: Bayesian 
demultiplexing of pooled single-cell RNA-seq data without genotype reference 
<https://www.biorxiv.org/content/10.1101/598748v1>`_. 
\ **bioRxiv** \ (2019): 598748.

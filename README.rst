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


Reference
=========

Yuanhua Huang, Davis J. McCarthy, and Oliver Stegle. `Vireo: Bayesian 
demultiplexing of pooled single-cell RNA-seq data without genotype reference 
<https://www.biorxiv.org/content/10.1101/598748v1>`_. 
\ **bioRxiv** \ (2019): 598748.

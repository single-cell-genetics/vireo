|PyPI| |Docs| |Build Status|

.. |PyPI| image:: https://img.shields.io/pypi/v/vireoSNP.svg
    :target: https://pypi.org/project/vireoSNP
.. |Docs| image:: https://readthedocs.org/projects/vireosnp/badge/?version=latest
   :target: https://vireoSNP.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/huangyh09/vireo.svg?branch=master
   :target: https://travis-ci.org/huangyh09/vireo


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

Vireo is available through PyPI_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install vireoSNP

Alternatively, you can download or clone this repository and type 
``python setup.py install`` to install. In either case, add ``--user`` if you 
don't have the permission as a root or for your Python environment.

For more instructions, see the installation_ manual.

.. _PyPI: https://pypi.org/project/vireoSNP
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

The vireoSNP python package offers a set of utilities functions and an  
executable command line `vireo` for donor deconvolution in any of these four 
situations:

**Mode 1:** without any genotype: 

.. code-block:: bash

   vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR

**Mode 2:** with genotype for all samples (specify tag ``-t``: GT, GP, or PL)

.. code-block:: bash

   vireo -c $CELL_DATA -d $DONOR_GT_FILE -o $OUT_DIR

**Mode 3:** with genotype for part of the samples (``N`` is different from the 
sample number in ``$DONOR_GT_FILE``)

.. code-block:: bash

   vireo -c $CELL_DATA -d $DONOR_GT_FILE -o $OUT_DIR -N $n_donor 

**Mode 4:** with genotype but not confident

.. code-block:: bash

   vireo -c $CELL_DATA -d $DONOR_GT_FILE -o $OUT_DIR --forceLearnGT

In modes 3 and 4 are less common, the algorithm will run mode 1 first and then 
match the estimated donor genotype and the given values (even partial). The 
matched given genotype will replace the estiamted ones as prior in the second 
run.

Note, the cell data (``$CELL_DATA``) via ``-c`` can be any of the following two 
formats:

* standard VCF file (compressed or uncompressed) with variants by cells
* a cellSNP output folder containing VCF for variants info and sparse matrices 
  `AD` and `DP`

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

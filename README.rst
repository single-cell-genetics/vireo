|PyPI| |Docs| |Build Status| |DOI|

.. |PyPI| image:: https://img.shields.io/pypi/v/vireoSNP.svg
    :target: https://pypi.org/project/vireoSNP
.. |Docs| image:: https://readthedocs.org/projects/vireosnp/badge/?version=latest
   :target: https://vireoSNP.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/single-cell-genetics/vireo.svg?branch=master
   :target: https://travis-ci.org/single-cell-genetics/vireo
.. |DOI| image:: https://zenodo.org/badge/187803798.svg
   :target: https://zenodo.org/badge/latestdoi/187803798


======================================================
vireo: donor deconvolution for pooled single-cell data
======================================================

Vireo: Variational Inference for Reconstructing Ensemble Origin by expressed 
SNPs in multiplexed scRNA-seq data. 

The name vireo follows the theme from cardelino_ (for clone deconvolution), 
while the Python package name is vireoSNP_ to aviod name confilict on PyPI.

.. _cardelino: https://github.com/PMBio/cardelino
.. _vireoSNP: https://pypi.org/project/vireoSNP


News
====
* All release notes can be found in doc/release.rst_.
* Notebook_ for subclone reconstructing with mitochrondrial mutations

.. _release.rst: https://github.com/single-cell-genetics/vireo/blob/master/doc/release.rst
.. _Notebook: https://vireosnp.readthedocs.io/en/latest/vireoSNP_clones.html


Installation
============

Vireo is available through PyPI_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install -U vireoSNP

Alternatively, you can install from this GitHub repository for latest (often 
development) version by following command line

.. code-block:: bash

  pip install -U git+https://github.com/single-cell-genetics/vireo

In either case, add ``--user`` if you don't have the write permission for your 
Python environment.

For more instructions, see the installation_ manual.

.. _PyPI: https://pypi.org/project/vireoSNP
.. _installation: https://vireoSNP.readthedocs.io/en/latest/install.html


Manual and examples
===================

The full manual is at https://vireoSNP.readthedocs.io 
It includes more details on installation, `demultiplex usage`_, and preprocess 
with genotyping_ cells.

Test example data is included in this repo and demos can be found in examples/demo.sh_.

Also, type ``vireo -h`` for all arguments with the version you are using.

.. _demultiplex usage: https://vireoSNP.readthedocs.io/en/latest/manual.html
.. _demo.sh: https://github.com/huangyh09/vireo/blob/master/examples/demo.sh
.. _genotyping: https://vireoSNP.readthedocs.io/en/latest/genotype.html


Reference
=========

Yuanhua Huang, Davis J. McCarthy, and Oliver Stegle. `Vireo: Bayesian 
demultiplexing of pooled single-cell RNA-seq data without genotype reference 
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2>`_. 
\ **Genome Biology** \ 20, 273 (2019)

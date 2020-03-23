|PyPI| |Docs| |Build Status| |DOI|

.. |PyPI| image:: https://img.shields.io/pypi/v/vireoSNP.svg
    :target: https://pypi.org/project/vireoSNP
.. |Docs| image:: https://readthedocs.org/projects/vireosnp/badge/?version=latest
   :target: https://vireoSNP.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/PMBio/vireo.svg?branch=master
   :target: https://travis-ci.org/PMBio/vireo
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
* vireoSNP_donors.ipynb_ giving example on donor deconvolution manually
* vireoSNP_clones.ipynb_ giving example clone reconstruction on mitochondral 
  mutations

.. _release.rst: https://github.com/single-cell-genetics/vireo/blob/master/doc/release.rst
.. _vireoSNP_donors.ipynb: https://github.com/single-cell-genetics/vireo/blob/master/examples/vireoSNP_donors.ipynb
.. _vireoSNP_clones.ipynb: https://github.com/single-cell-genetics/vireo/blob/master/examples/vireoSNP_clones.ipynb


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

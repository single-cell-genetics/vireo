|PyPI| |Docs| |Build Status| |DOI|

.. |PyPI| image:: https://img.shields.io/pypi/v/vireoSNP.svg
    :target: https://pypi.org/project/vireoSNP
.. |Docs| image:: https://readthedocs.org/projects/vireosnp/badge/?version=latest
   :target: https://vireoSNP.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/single-cell-genetics/vireo.svg?branch=master
   :target: https://travis-ci.org/single-cell-genetics/vireo
.. |DOI| image:: https://zenodo.org/badge/187803798.svg
   :target: https://zenodo.org/badge/latestdoi/187803798

====
Home
====

.. :Author: Yuanhua Huang
.. :Version: 0.2.0
.. :Last viewed: Jun 30, 2019

About Vireo
===========

This documentation gives an introduction and usage manual of Vireo (Variational 
inference for reconstructing ensemble origins), a Bayesian method to demultiplex
pooled scRNA-seq data with or without genotype reference.

Vireo is primarily designed for demultiplexing cells into donors by modelling of
expressed alleles. It supports a variety of settings of donor genotype (from
entirely missing, to partially missing, to fully observed). See more details in
`manual`_ section.

As a general cell clustering methods by allelic ratio (equivalent to genotyping),
Vireo is applicable for more settings besides donor demultiplexing, including
reconstruction of somatic clones, see `vireoSNP_clones.ipynb`_ for example on 
mitochondral mutations.

.. _manual: https://vireosnp.readthedocs.io/en/latest/manual.html


Notebooks for interactive analysis
==================================
Here are some notebooks for interactive analysis. Usually, you only need to use 
the command line to perform donor deconvolution, but you may refer to some of 
these notebooks for additional analysis.

**donors**: `vireoSNP_donors.ipynb`_ gives example on donor deconvolution 
manually. The `vireo` command line does this job automatically.

**donors**: `donor_match.ipynb`_ gives example on aligning donors to other
omics data or other batches

**clones**: `vireoSNP_clones.ipynb`_ gives example on clone reconstruction on 
mitochondral mutations

.. _donor_match.ipynb: https://nbviewer.jupyter.org/github/single-cell-genetics/vireo/blob/master/examples/donor_match.ipynb
.. _vireoSNP_donors.ipynb: https://nbviewer.jupyter.org/github/single-cell-genetics/vireo/blob/master/examples/vireoSNP_donors.ipynb
.. _vireoSNP_clones.ipynb: https://nbviewer.jupyter.org/github/single-cell-genetics/vireo/blob/master/examples/vireoSNP_clones.ipynb



Quick Resources
===============

**Latest version on GitHub**
https://github.com/single-cell-genetics/vireo

**Scripts for simulation**
https://github.com/single-cell-genetics/vireo/tree/master/simulate

**All releases**
https://pypi.org/project/vireoSNP/#history


Issue reports
=============
If you find any error or suspicious bug, we will appreciate your report.
Please write them in the github issues: 
https://github.com/single-cell-genetics/vireo/issues


References
==========

Yuanhua Huang, Davis J. McCarthy, and Oliver Stegle. `Vireo: Bayesian 
demultiplexing of pooled single-cell RNA-seq data without genotype reference 
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2>`_. 
\ **Genome Biology** \ 20, 273 (2019)


.. toctree::
    :caption: Main
    :maxdepth: 1
    :hidden:
    
    install
    manual
    genotype
    API
    release


.. toctree::
    :caption: Notebooks
    :maxdepth: 1
    :hidden:
 
    vireoSNP_clones

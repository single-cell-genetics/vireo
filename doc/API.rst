API
===

.. automodule:: vireoSNP

Import vireoSNP as::

   import vireoSNP


Commands
--------

* ``vireo``: see manual_
* ``GTbarcode``: see manual_

.. _manual: https://vireosnp.readthedocs.io/en/latest/manual.html

Read / Load
-----------

.. autofunction:: vireoSNP.read_cellSNP

.. autofunction:: vireoSNP.read_vartrix


VCF processing
--------------

Load VCF to matrices
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: vireoSNP.vcf.load_VCF

Parse genotype probablity to tenseor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: vireoSNP.vcf.parse_donor_GPb

.. autofunction:: vireoSNP.vcf.match_VCF_samples


Plotting
--------

Heatmap plot
~~~~~~~~~~~~

.. autofunction:: vireoSNP.plot.heat_matrix
   

Annotated heatmap plot
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: vireoSNP.plot.anno_heat



Vireo Object
------------

Objects of type :class:`~vireoSNP.Vireo` allow clustering cells by allelic ratio

.. autoclass:: vireoSNP.Vireo
   :members: __init__, set_initial, set_prior, fit


BinomMixtureVB Object
---------------------

Objects of type :class:`~vireoSNP.BinomMixtureVB` for clustering with binomial
mixture model

.. autoclass:: vireoSNP.BinomMixtureVB
   :members: __init__, set_initial, set_prior, fit


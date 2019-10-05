=======
History
=======

Release v0.1.6 (05/10/2019)
===========================
* Fix a bug when variants in donor genotype are not in cell vcf file

Release v0.1.5 (28/09/2019)
===========================
* Support genotype barcode generation

Release v0.1.4 (22/09/2019)
===========================
* Support that the case that input GT is larger than wanted `n_donor` 
* Clarify the structure in vireo_flock: 1) warm-up for multiple initials or 
  extra donors; 2) pre-step to subset or fill up the genotype prior; 3) the main
  run.
* Provide more options in the warm-up step to search donors from extra clusters.
  Before, it only uses the size of the donor. Now, the genotype distance can be
  used to search the K donors with furthest genotype distance.

Release v0.1.3 (30/08/2019)
===========================
* Support vartrix sparse matrices as input
* Change --amplifyK to --extraDonor for extra donors in initial search
* Fixed the bug for --noDoublet
* Fixed a bug for unassigned
* Minor update of figure output
* Updated the submoduals for easier import

Release v0.1.2 (15/07/2019)
===========================
* Support sparse matrices as input (for cellSNP directory with `-O`)
* Plot the distance between genotype probability between estimated samples
* Upgrade the manual, including the usage of simulation (readme in the 
  simulation folder of GitHub repo)

Release v0.1.1 (30/06/2019)
===========================
* A completed version for all planned features
* Donor deconvolution with supporting multiple modes:
  1) without genotype
  2) with genotype for all samples
  3) with genotype for part of the samples
  4) with genotype but not confident
* Manual for installation, usage, and preprocessing
* Release test data sets
* vireoSNP is available on PyPI, try it `pip install vireoSNP`

Release v0.1.0 (24/06/2019)
===========================
* reimplementation of vireo in Python (orignal in cardelino R package)
* Initial release with limited features

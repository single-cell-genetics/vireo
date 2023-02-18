=======
History
=======

Development on GitHub
=====================


Release v0.5.8 (18/02/2023)
===========================
* fix issue with None in match() and match_SNPs() when w/ chr w/o chr have partial match
* minor optimise the codes for vireoSNP.utils.vcf_utils.parse_donor_GPb
* add a snp_gene_match() function
* fix issue in write_VCF() when there is no samples

Release v0.5.7 (24/03/2022)
===========================
* fix the issue when output_dir is not given
* add doc for read VCF files

Release v0.5.6 (07/04/2021)
===========================
* fix a bug for detecting unsupported genotyep tag
* add a wrap function to compare samples in two VCF files
* add doublet_logLikRatio in `donor_ids.tsv` for extra indicators of doublets
* update documentation with supporting notebook `vireoSNP_clones.ipynb`
* update API

Release v0.5.5 (28/03/2021)
===========================
* update notebook `vireoSNP_clones.ipynb`
* update API

Release v0.5.4 (28/03/2021)
===========================
* introduce log likelihood ratio for detecting ambient RNAs
* support ambient RNAs from a mixture of all donors or only two donors
* introduce multiple processes for multiple initializations
* introduce ELBO_gain for selecting variants
* For donor_ids.tsv, the doublet_prob change from sum to max

Release v0.5.3 (28/03/2021)
===========================
* support detection of ambient RNAs, alternative way for doublet detection

Release v0.5.0 (09/02/2021)
===========================
* support support numpy.ndarray and automatically change to sparse matrix
* fix the sign of KL
* fix a minor bug on --noDoublet setting
* add --cellRange to subset the input cells for less memory
* update anno_heat() plotting
* add get_confusion() for results comparison and plotting

Release v0.4.2 (14/11/2020)
===========================
* fix the donor names when N < donors_in_GT
* change the suggestion from cellSNP python to C version (cellsnp-lite)

Release v0.4.1 (18/05/2020)
===========================
* add likelihood ratio test the differential donor abundance in bulk RNA-seq
  data
* set the interpolation='none' for plt.imshow

Release v0.4.0 (19/04/2020)
===========================
* add vireoBulk for demultiplexing in bulk RNA-seq data

Release v0.3.2 (10/04/2020)
===========================
* support donor variant match between with and without "chr" prefix

Release v0.3.1 (25/03/2020)
===========================
* replace greed_match to optimal_match for aligning donors via genotype

Release v0.3.0 (23/03/2020)
===========================
* Rewrite the Vireo in the object-oriented way for easier upgrading and adding 
  new features
* Now support fix the dispersion of the theta posterior distribution
* Change `delay_fit_theta` as an augument. It's often useful for donor
  deconvolution, but may not ideal for clonal inference, where theta can be very 
  different from our expectation due to ASE or copy numbers

Release v0.2.3 (22/03/2020)
===========================
* Fix a minor bug in donor_select()

Release v0.2.2 (21/03/2020)
===========================
* Change GP_prob's shape from (n_var, n_GT, n_donor) to (n_var, n_donor, n_GT)
* Restructure the codes for further upgrading
* Minor fix the GT_plot xlim and ylim

Release v0.2.1 (30/01/2020)
===========================
* Fix a bug when the donors in the input GT is smaller than donors in the pooled
  scRNA-seq. The sample id is now corrected.

Release v0.2.0 (28/01/2020)
===========================
* Support SNP specific allelic ratio, namely theta parameters. Note, SNP based 
  ASE mode requires a much stronger prior on theta to avoid overfitting, as each
  variant has very low number of reads. 
* Change the default extraDonor to 0.
* Provide examples/vireoSNP_usage.ipynb for using vireoSNP as a Python module 
  for general cell clustering based on allelic ratio.

Release v0.1.8 (29/10/2019)
===========================
* Further fix the bug when variants in donor genotype are not in cell vcf file

Release v0.1.7 (05/10/2019)
===========================
* Support donor genotype vcf file with different FORMAT for different variants

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

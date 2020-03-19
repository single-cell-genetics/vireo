=============================
Synthetic mixture of 10x data
=============================

In this folder, we provide the script to generate synthetic mixture of 10x data
from multiple individuals.

The python script is `synth_pool.py`_, which can mix multiple BAM files from 
10x cellranger outputs. There are a few features supported for the synthetic 
mixture:

* Pooling multiple BAM files (each BAM file represents an individual) by using
  ``--samFiles`` for comma separated BAM files
* For each BAM file, only keeping the reads for given cells by using 
  ``--barcodeFiles`` for comma separated cell barcode files. This should be the 
  same number of files to the BAM files. The sample id will be added as suffix 
  to the cell barcodes, e.g., xxx-3 for the 3rd sample. If it's a doublet it 
  will be xxx-3D if its a doublet from two samples or xxx-3S if it's from the 
  same sample 3.
* For all BAM files, only keep the reads covering the given variants by using 
  ``--regionFile`` in VCF format
* Extra cells are added to the same number of singletons to generate doublets. 
  The number of doublets are defined by the total cells and the doublet rate
  ``--doubletRate``
* The number of cells for each sample ``--nCELL``, cell counts are equally 
  distributed across all samples, but can make one of the sample as a minor 
  sample with using ``--minorSAMPLE`` for its ratio to other samples.

.. _synth_pool.py: https://github.com/huangyh09/vireo/blob/master/simulate/synth_pool.py

For details of all arguments: type ``python synth_pool.py  -h``

.. code-block:: html

   Usage: synth_pool.py [options]

   Options:
    -h, --help            show this help message and exit
    -s SAM_FILES, --samFiles=SAM_FILES
                            Input sam files, comma separated.
    -b BARCODES_FILES, --barcodeFiles=BARCODES_FILES
                            Input barcode files, comma separated.
    -r REGION_FILE, --regionFile=REGION_FILE
                            Input SNP list.
    -d DOUBLET_RATE, --doubletRate=DOUBLET_RATE
                            Doublet rate [default: n/100000]
    -o OUT_DIR, --outDir=OUT_DIR
                            Directory for output files: pooled.bam and
                            barcodes_pool.tsv.
    -p NPROC, --nproc=NPROC
                            Number of subprocesses [default: 1]

    Cell barcodes sampling:
        --nCELL=N_CELL      The number of cells in each sample [default: none]
        --minorSAMPLE=MINOR_SAMPLE
                            Ratio size of minor sample [default: 1.0]
        --randomSEED=RANDOM_SEED
                            The random seed in numpy [default: none]

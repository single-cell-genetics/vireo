# vireoSNP - donor deconvolution for multiplexed scRNA-seq data
# Author: Yuanhua Huang
# Date: 24-06-2019

import os
import sys
import time
import subprocess
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from .version import __version__
from .utils.vireo_model import vireo_core
from .utils.io_utils import write_donor_id
from .utils.vcf_utils import load_VCF, write_VCF
from .utils.vcf_utils import read_sparse_GeneINFO, GenoINFO_maker


START_TIME = time.time()

def show_progress(RV=None):
    return RV

def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--cellFile", "-c", dest="cell_file", default=None,
        help=("The cell genotyping file in VCF or H5 format."))
    parser.add_option("--nDonor", "-N", type="int", dest="n_donor", 
        default=None, help="Number of donors to infer [default: %default]")
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Dirtectory for output files [default: $cellFilePath/vireo]"))

    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--donorFile", "-d", dest="donor_file", default=None,
        help=("The donor genotyping file in VCF or H5 format [default: NA]"))
    group1.add_option("--GTtag", "-t", dest="GT_tag", default='PL',
        help=("The tag for donor genotype: GT, GP, PL [default: %default]"))
    group1.add_option("--noDoublet", dest="no_doublet", action="store_true", 
        default=False, help="If use, not checking doublets.")
    group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    group1.add_option("--randSeed", type="int", dest="rand_seed", default=None,
        help="Seed for random initialization [default: %default]")
    
    parser.add_option_group(group1)
    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        print("Welcome to vireoSNP v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)

    ## out directory
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $cellFilePath/vireo.")
        out_dir = os.path.dirname(os.path.abspath(options.cell_file)) + "/vireo"
    elif os.path.dirname(options.out_dir) == "":
        out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ## input data
    n_donor = options.n_donor
    donor_names = ['donor%d' %x for x in range(n_donor)]
        
    if options.cell_file is None:
        print("Error: need cell file in vcf or hdf5 format.")
        sys.exit(1)
    else:
        print("[vireo] Loading cell VCF file ...")
        cell_vcf = load_VCF(options.cell_file)
        cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
        n_vars = np.array(np.sum(cell_dat['DP'] > 0, axis=0)).reshape(-1)
        print("[vireo] Loading cell VCF file, done.")

    ## run vireo model (try multiple initializations)
    res_vireo = vireo_core(cell_dat['AD'], cell_dat['DP'], n_donor=n_donor, 
        random_seed=options.rand_seed)
    print("[vireo] VB lower bound: %.2f" %res_vireo['LB_list'][-1])  
    
    ## save donor id for each cell
    write_donor_id(out_dir, donor_names, cell_vcf['samples'], n_vars,
        res_vireo['ID_prob'], res_vireo['doublet_prob'])

    # ## save inferred donor genotype
    # donor_vcf_out = cell_vcf
    # donor_vcf_out['samples'] = donor_names
    # donor_vcf_out['GenoINFO'] = GenoINFO_maker(res_vireo['GT_prob'], 
    #     res_vireo['AD_donor'], res_vireo['DP_donor'])
    # write_VCF(out_dir + "/GT_donors.vireo.vcf.gz", donor_vcf_out)
    
    run_time = time.time() - START_TIME
    print("[vireo] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                   run_time % 60))
    
        
if __name__ == "__main__":
    main()

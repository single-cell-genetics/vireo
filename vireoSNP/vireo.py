# vireoSNP - donor deconvolution for multiplexed scRNA-seq data
# Author: Yuanhua Huang
# Date: 24-06-2019

import os
import sys
import time
import subprocess
import numpy as np
import multiprocessing
from scipy.io import mmread
from optparse import OptionParser, OptionGroup

from .version import __version__
from .utils.vireo_base import match
from .utils.vireo_model import vireo_core, vireo_flock
from .utils.io_utils import write_donor_id, plot_GT
from .utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb
from .utils.vcf_utils import read_sparse_GeneINFO, GenoINFO_maker


START_TIME = time.time()

def show_progress(RV=None):
    return RV

def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--cellData", "-c", dest="cell_data", default=None,
        help=("The cell genotype file in VCF format or cellSNP folder with "
              "sparse matrices."))
    parser.add_option("--nDonor", "-N", type="int", dest="n_donor", 
        default=None, help=("Number of donors to demultiplex; can be larger "
        "than provided in donor_file"))
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Dirtectory for output files [default: $cellFilePath/vireo]"))

    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--noDoublet", dest="no_doublet", action="store_true", 
        default=False, help="If use, not checking doublets.")
    group1.add_option("--donorFile", "-d", dest="donor_file", default=None,
        help=("The donor genotype file in VCF format. Please filter the sample "
        "and region with bcftools -s and -R first!"))
    group1.add_option("--genoTag", "-t", dest="geno_tag", default='PL',
        help=("The tag for donor genotype: GT, GP, PL [default: %default]"))
    group1.add_option("--nInit", "-M", type="int", dest="n_init", default=None,
        help="Number of random initializations [default: 2 (GT) or 50 (no GT)]")
    group1.add_option("--amplifyK", type=float, dest="K_amplify", default=None,
        help="Pre-cluster with amplified K [default: 1.0 (GT) or 1.2 (no GT)]")
    group1.add_option("--forceLearnGT", dest="force_learnGT", default=False, 
        action="store_true", help="If use, treat donor GT as prior only.")
    group1.add_option("--noPlot", dest="no_plot", default=False, 
        action="store_true", help="If use, turn off plotting GT distance.")
    group1.add_option("--randSeed", type="int", dest="rand_seed", default=None,
        help="Seed for random initialization [default: %default]")
    # group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
    #     help="Number of subprocesses [default: %default]")
    
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

    ## input data (VCF.gz or a folder with sparse matrices)
    if options.cell_data is None:
        print("Error: need cell file in vcf or hdf5 format.")
        sys.exit(1)
    elif os.path.isdir(os.path.abspath(options.cell_data)):
        print("[vireo] Loading cell folder ...")
        cell_vcf = load_VCF(options.cell_data + "/cellSNP.base.vcf.gz", 
            load_sample=False)
        cell_vcf['samples'] = np.genfromtxt(options.cell_data + 
            "/cellSNP.samples.tsv", dtype=str)

        cell_dat = {}
        cell_dat['AD'] = mmread(options.cell_data + "/cellSNP.tag.AD.mtx").tocsc()
        cell_dat['DP'] = mmread(options.cell_data + "/cellSNP.tag.DP.mtx").tocsc()
        print("[vireo] Loading cell foder, done.")
    else:
        print("[vireo] Loading cell VCF file ...")
        cell_vcf = load_VCF(options.cell_data)
        cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
        print("[vireo] Loading cell VCF file, done.")
    n_vars = np.array(np.sum(cell_dat['DP'] > 0, axis=0)).reshape(-1)

    ## input donor genotype
    n_donor = options.n_donor
    if options.donor_file is not None:
        print("[vireo] Loading donor VCF file ...")
        donor_vcf = load_VCF(options.donor_file, sparse=False)
        if (options.geno_tag not in donor_vcf['GenoINFO']):
            print("[vireo] No " + options.geno_tag + " tag in donor genotype; " 
                "please try another tag for genotype, e.g., GT")
            print("        %s" %options.donor_file)
            sys.exit(1)
        donor_GPb = parse_donor_GPb(donor_vcf['GenoINFO'][options.geno_tag], 
            options.geno_tag)
        
        mm_idx = match(cell_vcf['variants'], donor_vcf['variants'])
        idx1 = np.where(mm_idx != None)[0]
        idx2 = mm_idx[idx1].astype(int)

        cell_dat['AD'] = cell_dat['AD'][idx1, :]
        cell_dat['DP'] = cell_dat['DP'][idx1, :]
        donor_GPb = donor_GPb[idx2, :, :]
        print("[vireo] Loading donor VCF file, done.")

        if n_donor is None or n_donor <= donor_GPb.shape[2]:
            n_donor = donor_GPb.shape[2]
            donor_names = donor_vcf['samples']
            learn_GT = False
        else:
            learn_GT = True
            donor_names = (donor_vcf['samples'] + 
                ['donor%d' %x for x in range(donor_GPb.shape[2], n_donor)])
    else:
        learn_GT = True
        donor_GPb = None
        donor_names = ['donor%d' %x for x in range(n_donor)]

    if options.force_learnGT:
        learn_GT = True
    if options.n_init is None:
        n_init = 50 if learn_GT else 2
    else:
        n_init = options.n_init
    if options.K_amplify is None:
        K_amplify = 1.2 if learn_GT else 1.0
    else:
        K_amplify = options.K_amplify

    ## run vireo model (try multiple initializations)
    print("[vireo] Demultiplex %d cells to %d donors with %d variants." %(
        cell_dat['AD'].shape[1], n_donor, cell_dat['AD'].shape[0]))
    res_vireo = vireo_flock(cell_dat['AD'], cell_dat['DP'], n_donor=n_donor, 
        GT_prior=donor_GPb, learn_GT=learn_GT, n_init=n_init, 
        K_amplify=K_amplify, random_seed=options.rand_seed)

    ## save donor id for each cell
    write_donor_id(out_dir, donor_names, cell_vcf['samples'], n_vars,res_vireo)

    if options.no_plot == False:
        idx = np.array(np.sum(cell_dat['DP'], axis=1) > (3*n_donor)).reshape(-1)
        if learn_GT and donor_GPb is not None:
            plot_GT(out_dir, res_vireo['GT_prob'][idx, :, :], donor_names, 
                    donor_GPb[idx, :, :], donor_vcf['samples'])
        else:
            plot_GT(out_dir, res_vireo['GT_prob'][idx, :, :], donor_names)

    # ## save inferred donor genotype
    if learn_GT:
        donor_vcf_out = cell_vcf
        donor_vcf_out['samples'] = donor_names
        donor_vcf_out['GenoINFO'] = GenoINFO_maker(res_vireo['GT_prob'], 
            cell_dat['AD'] * res_vireo['ID_prob'], 
            cell_dat['DP'] * res_vireo['ID_prob'])
        write_VCF(out_dir + "/GT_donors.vireo.vcf.gz", donor_vcf_out)
    
    run_time = time.time() - START_TIME
    print("[vireo] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                    run_time % 60))
    
        
if __name__ == "__main__":
    main()

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
from .utils.vireo_base import match, greed_match
from .utils.vireo_wrap import vireo_wrap

from .plot.base_plot import plot_GT
from .utils.io_utils import write_donor_id, read_cellSNP, read_vartrix
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

    group0 = OptionGroup(parser, "Optional input files")
    group0.add_option("--vartrixData", dest="vartrix_data", default=None,
        help=("The cell genotype files in vartrix outputs (three/four files, "
              "comma separated): alt.mtx,ref.mtx,barcodes.tsv,SNPs.vcf.gz. "
              "This will suppress cellData argument."))
    group0.add_option("--donorFile", "-d", dest="donor_file", default=None,
        help=("The donor genotype file in VCF format. Please filter the sample "
        "and region with bcftools -s and -R first!"))
    group0.add_option("--genoTag", "-t", dest="geno_tag", default='PL',
        help=("The tag for donor genotype: GT, GP, PL [default: %default]"))
    
    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--noDoublet", dest="no_doublet", action="store_true", 
        default=False, help="If use, not checking doublets.")
    group1.add_option("--nInit", "-M", type="int", dest="n_init", default=50,
        help=("Number of random initializations, when GT needs to learn "
        "[default: %default]"))
    group1.add_option("--extraDonor", type=int, dest="n_extra_donor", 
        default=0, help=("Number of extra donor in pre-cluster, when GT "
        "needs to learn [default: %default]"))
    group1.add_option("--extraDonorMode", dest="extra_donor_mode", 
        default="distance", help=("Method for searching from extra donors. "
        "size: n_cell per donor; distance: GT distance between donors "
        "[default: %default]"))
    group1.add_option("--forceLearnGT", dest="force_learnGT", default=False, 
        action="store_true", help="If use, treat donor GT as prior only.")
    group1.add_option("--ASEmode", dest="ASE_mode", default=False, 
        action="store_true", help="If use, turn on SNP specific allelic ratio.")
    group1.add_option("--noPlot", dest="no_plot", default=False, 
        action="store_true", help="If use, turn off plotting GT distance.")
    group1.add_option("--randSeed", type="int", dest="rand_seed", default=None,
        help="Seed for random initialization [default: %default]")
    # group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
    #     help="Number of subprocesses [default: %default]")
    
    parser.add_option_group(group0)
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
    if options.cell_data is None and options.vartrix_data is None:
        print("Error: need cell data in vcf file, or cellSNP output folder, or "
              "vartrix's alt.mtx,ref.mtx,barcodes.tsv.")
        sys.exit(1)
    elif options.vartrix_data is not None:
        print("[vireo] Loading vartrix files ...")
        vartrix_files = options.vartrix_data.split(",")
        if len(vartrix_files) < 3 or len(vartrix_files) > 4:
            print("Error: vartrixData requires 3 or 4 comma separated files")
            sys.exit(1)
        elif len(vartrix_files) == 3:
            vartrix_files.append(None)

        cell_dat = read_vartrix(vartrix_files[0], vartrix_files[1], 
                                vartrix_files[2], vartrix_files[3])
    elif os.path.isdir(os.path.abspath(options.cell_data)):
        print("[vireo] Loading cell folder ...")
        cell_dat = read_cellSNP(options.cell_data)
    else:
        print("[vireo] Loading cell VCF file ...")
        cell_vcf = load_VCF(options.cell_data, biallelic_only=True)
        cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
        for _key in ['samples', 'variants', 'FixedINFO', 'contigs', 'comments']:
            cell_dat[_key] = cell_vcf[_key]

    ## input donor genotype
    n_donor = options.n_donor
    if options.donor_file is not None:
        if "variants" not in cell_dat.keys():
            print("No variants information is loaded, please provide base.vcf.gz")
            sys.exit(1)
        
        print("[vireo] Loading donor VCF file ...")
        donor_vcf = load_VCF(options.donor_file, biallelic_only=True, 
                             sparse=False, format_list=[options.geno_tag])
        if (options.geno_tag not in donor_vcf['GenoINFO']):
            print("[vireo] No " + options.geno_tag + " tag in donor genotype; " 
                "please try another tag for genotype, e.g., GT")
            print("        %s" %options.donor_file)
            sys.exit(1)
        donor_GPb = parse_donor_GPb(donor_vcf['GenoINFO'][options.geno_tag], 
            options.geno_tag)
        
        mm_idx = match(cell_dat['variants'], donor_vcf['variants'])
        mm_idx = mm_idx.astype(float)
        idx1 = np.where(mm_idx == mm_idx)[0] #remove None
        # TODO: check when chr is not compatible! given warning.
        if len(idx1) == 0:
            print("[vireo] warning: no variants matched to donor VCF, " + 
                  "please check chr format!")
        else:
            print("[vireo] %d out %d variants matched to donor VCF" 
                  %(len(idx1), len(cell_dat['variants'])))
        idx2 = mm_idx[idx1].astype(int)

        donor_GPb = donor_GPb[idx2, :, :]
        cell_dat['AD'] = cell_dat['AD'][idx1, :]
        cell_dat['DP'] = cell_dat['DP'][idx1, :]
        cell_dat["variants"]  = [cell_dat["variants"][x] for x in idx1]
        for _key in cell_dat["FixedINFO"].keys():
            cell_dat["FixedINFO"][_key] = [
                cell_dat["FixedINFO"][_key][x] for x in idx1]

        if n_donor is None or n_donor == donor_GPb.shape[1]:
            n_donor = donor_GPb.shape[1]
            donor_names = donor_vcf['samples']
            learn_GT = False
        elif n_donor < donor_GPb.shape[1]:
            learn_GT = False
            donor_names = ['donor%d' %x for x in range(n_donor)]
        else:
            learn_GT = True
            donor_names = (donor_vcf['samples'] + 
                ['donor%d' %x for x in range(donor_GPb.shape[1], n_donor)])
    else:
        learn_GT = True
        donor_GPb = None
        donor_names = ['donor%d' %x for x in range(n_donor)]

    n_vars = np.array(np.sum(cell_dat['DP'] > 0, axis=0)).reshape(-1)

    if options.force_learnGT:
        learn_GT = True
    
    # extra donor for initial search, only for learn_GT
    n_extra_donor = 0
    if learn_GT:
        if options.n_extra_donor is None or options.n_extra_donor == "None":
            n_extra_donor = int(round(np.sqrt(n_donor)))
        else:
            n_extra_donor = options.n_extra_donor        
    
    # number of initials, only for learn_GT
    n_init = options.n_init if learn_GT else 1
    
    check_doublet = options.no_doublet == False

    ## run vireo model (try multiple initializations)
    print("[vireo] Demultiplex %d cells to %d donors with %d variants." %(
        cell_dat['AD'].shape[1], n_donor, cell_dat['AD'].shape[0]))
    res_vireo = vireo_wrap(cell_dat['AD'], cell_dat['DP'], n_donor=n_donor, 
        GT_prior=donor_GPb, learn_GT=learn_GT, n_init=n_init, 
        n_extra_donor=n_extra_donor, extra_donor_mode=options.extra_donor_mode,
        check_doublet=check_doublet, random_seed=options.rand_seed,
        ASE_mode=options.ASE_mode)

    if (n_donor is not None and 
        donor_GPb is not None and n_donor < donor_GPb.shape[2]):
        idx = greed_match(res_vireo['GT_prob'], donor_GPb)
        donor_names = [donor_vcf['samples'][x] for x in idx]

    ## save donor id for each cell
    write_donor_id(out_dir, donor_names, cell_dat['samples'], n_vars, res_vireo)

    if options.no_plot == False and options.vartrix_data is None:
        idx = np.array(np.sum(cell_dat['DP'], axis=1) > (3*n_donor)).reshape(-1)
        if learn_GT and donor_GPb is not None:
            plot_GT(out_dir, res_vireo['GT_prob'][idx, :, :], donor_names, 
                    donor_GPb[idx, :, :], donor_vcf['samples'])
        else:
            plot_GT(out_dir, res_vireo['GT_prob'][idx, :, :], donor_names)

    # ## save inferred donor genotype
    if learn_GT and 'variants' in cell_dat.keys():
        donor_vcf_out = cell_dat
        donor_vcf_out['samples'] = donor_names
        donor_vcf_out['GenoINFO'] = GenoINFO_maker(res_vireo['GT_prob'], 
            cell_dat['AD'] * res_vireo['ID_prob'], 
            cell_dat['DP'] * res_vireo['ID_prob'])
        write_VCF(out_dir + "/GT_donors.vireo.vcf.gz", donor_vcf_out)
    
    run_time = time.time() - START_TIME
    print("[vireo] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                    run_time % 60))
    print()
    
        
if __name__ == "__main__":
    main()

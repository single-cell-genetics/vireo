# CLI for mitoMut

import os
import sys
import time
import numpy as np
from scipy.io import mmread
from optparse import OptionParser, OptionGroup

from vireoSNP import __version__
from vireoSNP.utils.io_utils import read_cellSNP, read_vartrix, read_sparse_GeneINFO
from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb

from .mitoMut import MitoMut

START_TIME = time.time()

def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--cellData", "-c", dest="cell_data", default=None,
        help=("The cellSNP folder with AD and DP sparse matrices."))
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Dirtectory for output files [default: $cellData/mitoMut]"))
    
    group0 = OptionGroup(parser, "Alternative input formats")
    group0.add_option("--mtxData", "-m", dest="mtx_data", default=None,
        help=("The two mtx files for AD and DP matrices, comma separated"))
    parser.add_option("--vcfData", dest="vcf_data", default=None,
        help=("The cell genotype file in VCF format"))
    
    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--randSeed", type="int", dest="rand_seed", default=None,
        help="Seed for random initialization [default: %default]")
    group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    
    parser.add_option_group(group0)
    parser.add_option_group(group1)
    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        print("Welcome to mitoMut v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)

    ## out directory
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $cellFilePath/mitoMut/")
        out_dir = os.path.dirname(os.path.abspath(options.cell_file)) + "/mitoMut/"
    elif os.path.dirname(options.out_dir) == "":
        out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ## input data (VCF.gz or a folder with sparse matrices)
    if (options.cell_data is None and options.mtx_data is None and 
        options.vcf_data is None):
        print("Error: need cell data in vcf file, or cellSNP output folder, or "
              "matrix files: AD.mtx,DP.mtx")
        sys.exit(1)
    elif options.mtx_data is not None:
        print("[mitoMut] Loading matrix files ...")
        matrix_files = options.mtx_data.split(",")
        if len(matrix_files) != 2:
            print("Error: mtxData requires 2 comma separated files")
            sys.exit(1)
        cell_dat = {}
        cell_dat['AD'] = mmread(matrix_files[0])
        cell_dat['DP'] = mmread(matrix_files[1])
        cell_dat['variants'] = ['SNP%d' %(x + 1) for x in 
                                range(cell_dat['AD'].shape[0])]
        
    elif options.vcf_data is not None:
        print("[mitoMut] Loading cell VCF file ...")
        cell_vcf = load_VCF(options.vcf_data, biallelic_only=True)
        cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
        for _key in ['samples', 'variants', 'FixedINFO', 'contigs', 'comments']:
            cell_dat[_key] = cell_vcf[_key]
        
    else:
        #os.path.isdir(os.path.abspath(options.cell_data)):
        print("[mitoMut] Loading cell folder ...")
        cell_dat = read_cellSNP(options.cell_data)
        
    ## More options
    nproc = options.nproc
    
    
    ## Main functions
    mdphd = MitoMut(AD = cell_dat['AD'], DP = cell_dat['DP'], 
                    variant_names = cell_dat['variants'])
    
    df = mdphd.fit_deltaBIC(out_dir = out_dir, nproc = nproc, beta_mode = False)
    print(df)
    final = mdphd.filter(by='deltaBIC', threshold = 500, out_dir = out_dir)
   
    
    run_time = time.time() - START_TIME
    print("[mitoMut] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                  run_time % 60))
    print()
    
        
if __name__ == "__main__":
    main()


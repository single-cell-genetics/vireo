# GTbarcode - generator of genotype barcode for discriminatory variants 
# Author: Yuanhua Huang
# Date: 28-09-2019

import os
import sys
import numpy as np
from optparse import OptionParser, OptionGroup

from .version import __version__
from .plot.base_plot import minicode_plot
from .utils.variant_select import variant_select
from .utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb


def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--vcfFile", "-i", dest="vcf_file", default=None,
        help="The VCF file for genotype of samples")
    parser.add_option("--outFile", "-o", dest="out_file", default=None,
        help="Output file [default: $vcfFile/GTbarcode.tsv]")

    group0 = OptionGroup(parser, "Optional arguments")
    group0.add_option("--genoTag", "-t", dest="geno_tag", default='GT',
        help=("The tag for donor genotype: GT, GP, PL [default: %default]"))
    group0.add_option("--noHomoAlt", dest="no_homo_alt", default=False, 
        action="store_true", help="Filter out variants with homozygous ALT.")
    group0.add_option("--noPlot", dest="no_plot", default=False, 
        action="store_true", help="Turn off the plot for the barcode.")
    group0.add_option("--figSize", dest="fig_size", default="4,2", 
        help="Size for the output figure, comma separated [default: %default].")
    group0.add_option("--figFormat", dest="fig_format", default="png", 
        help="Format of output figure: png or pdf [default: %default].")
    group0.add_option("--randSeed", type="int", dest="rand_seed", default=None,
        help=("Seed for random pick variants with same information gain "
              "[default: %default]"))
    
    parser.add_option_group(group0)
    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        print("Welcome to GT barcode generator; Vireo v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)

    ## input data vcf.gz
    if options.vcf_file is None:
        print("Error: need genotype data in vcf file.")
        sys.exit(1)
    else:
        vcf_file = options.vcf_file

    ## out directory
    if options.out_file is None:
        print("Warning: no outFile provided, we use $vcfFile/GTbarcode.tsv")
        out_file = os.path.dirname(os.path.abspath(vcf_file)) + "/GTbarcode.tsv"
    else:
        out_file = options.out_file
    if not os.path.exists(os.path.dirname(out_file)):
        os.mkdir(os.path.dirname(out_file))


    ## Load VCF data
    geno_tag = options.geno_tag
    donor_vcf = load_VCF(vcf_file, sparse=False, biallelic_only=True)
    donor_GPb = parse_donor_GPb(donor_vcf['GenoINFO'][geno_tag], geno_tag)

    var_ids = np.array(donor_vcf["variants"])
    GT_vals = np.argmax(donor_GPb, axis = 2)
    sample_ids = donor_vcf['samples']

    INFO = donor_vcf["FixedINFO"]["INFO"]
    AD, DP, OTH = [], [], []
    for i in range(len(INFO)):
        if INFO[i].count("AD=") == 0:
            AD.append(0)
        else: 
            AD.append(float(INFO[i].split("AD=")[1].split(";")[0]))
            
        if INFO[i].count("DP=") == 0:
            DP.append(0)
        else: 
            DP.append(float(INFO[i].split("DP=")[1].split(";")[0]))
            
        if INFO[i].count("OTH=") == 0:
            OTH.append(0)
        else: 
            OTH.append(float(INFO[i].split("OTH=")[1].split(";")[0]))
    AD, DP, OTH = np.array(AD), np.array(DP), np.array(OTH)

    ## filetering
    idx = (DP > 20) * (OTH / DP < 0.05)
    if options.no_homo_alt:
        idx *= np.max(GT_vals, axis=1) < 2

    AD, DP, OTH = AD[idx], DP[idx], OTH[idx]
    var_ids, GT_vals = var_ids[idx], GT_vals[idx, :]

    res_barcodes = variant_select(GT_vals, DP, rand_seed=options.rand_seed)
    fid = open(out_file, "w")
    #fid.writelines("\n".join(["barcodes"] + sample_ids) + "\n")
    fid.writelines("\t".join(["variants"] + sample_ids) + "\n")
    for i in res_barcodes[2]:
        line_list = [var_ids[i]] + ["%d" %x for x in GT_vals[i, :]]
        fid.writelines("\t".join(line_list) + "\n")
    fid.close()

    ## plot
    if options.no_plot == False:
        fig_size = np.array(options.fig_size.split(","), float)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(fig_size[0], fig_size[1]), dpi=300)
        minicode_plot(res_barcodes[1], var_ids[res_barcodes[2]], 
            donor_vcf['samples'])
        plt.tight_layout()
        fig.savefig(".".join(out_file.split(".")[:-1]) + "." + 
                    options.fig_format)


if __name__ == "__main__":
    main()

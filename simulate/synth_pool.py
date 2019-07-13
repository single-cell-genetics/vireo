# Sythetic mixture of bam files from multiple samples
# Author: Yuanhua Huang
# Date: 15-06-2019

import os
import sys
import pysam
import itertools
import numpy as np
import subprocess
import multiprocessing
from optparse import OptionParser, OptionGroup
from cellSNP.utils.vcf_utils import load_VCF
from cellSNP.utils.pileup_utils import check_pysam_chrom


def show_progress(RV=None):
    return RV


def sample_barcodes(barcodes, n_cell_each=1000, minor_sample=1.0, seed=None):
    """
    generate cell barcodes by down sampling
    """
    if seed is not None:
        np.random.seed(seed)
    for ss in range(len(barcodes)):
        if len(barcodes[ss]) < n_cell_each:
            print("Error in sample_barcodes: input sample has fewer cell "
                  "barcodes than n_cell_each.")
            sys.exit(1)
        barcodes[ss] = list(np.random.permutation(barcodes[ss])[:n_cell_each])
    barcodes[0] = barcodes[0][:round(minor_sample * n_cell_each)]
    return barcodes


def pool_barcodes(barcodes, out_dir, doublet_rate=None, sample_suffix=True, 
    seed=None):
    """
    Update cell barcodes with sample id and add doublets.
    Note, barcodes is a list of multiple samples, each  
    sample has a list of barcodes.
    """
    if seed is not None:
        np.random.seed(seed)
    if sample_suffix:
        barcodes_out = []
        for ss in range(len(barcodes)):
            barcodes_out.append([x[:-1]+str(ss+1) for x in barcodes[ss]])
    else:
        barcodes_out = barcodes.copy()
    barcodes_flat = list(itertools.chain(*barcodes_out))
            
    n_cells = len(barcodes_flat)
    if doublet_rate is None:
        doublet_rate = n_cells / 100000.0
    elif doublet_rate < 0 or doublet_rate > 1:
        print("Error: doublet rate needs to be between 0 and 1.")
        sys.exit(1)
    if doublet_rate == 0:
        n_doublets = 0
    else:
        n_doublets = round(n_cells / (1 + 1 / doublet_rate))
        
    print(n_cells, n_doublets)

    perm_idx = np.random.permutation(n_cells)
    for ii in range(n_doublets):
        if (barcodes_flat[perm_idx[ii]].split("-")[1] == 
            barcodes_flat[perm_idx[ii + n_doublets]].split("-")[1]):
            _barcode = barcodes_flat[perm_idx[ii]] + "S"
        else:
            _barcode = barcodes_flat[perm_idx[ii]] + "D"
        barcodes_flat[perm_idx[ii]] = _barcode
        barcodes_flat[perm_idx[ii + n_doublets]] = _barcode

    start_idx = 0
    for ss in range(len(barcodes_out)):
        _n_cell = len(barcodes_out[ss])
        barcodes_out[ss] = barcodes_flat[start_idx: start_idx + _n_cell]
        start_idx += _n_cell

    ## save new cell barcodes
    fid = open(out_dir + "/barcodes_pool.tsv", "w")
    for _barcode in np.unique(barcodes_flat):
        fid.writelines(_barcode + "\n")
    fid.close()

    fid = open(out_dir + "/cell_info.tsv", "w")
    fid.writelines("CB_pool\tCB_origin\tSample_id\n")
    for ss in range(len(barcodes_out)):
        for ii in range(len(barcodes_out[ss])):
            _out = [barcodes_out[ss][ii], barcodes[ss][ii], str(ss + 1)]
            fid.writelines("\t".join(_out) + "\n")
    fid.close()
    return barcodes_out


def fetch_reads(samFile_list, chroms, positions, outbam, 
                barcodes_in, barcodes_out=None, cell_tag='CB'):
    """
    """
    samFile_list = [check_pysam_chrom(x, chroms[0])[0] for x in samFile_list]
    outbam = pysam.AlignmentFile(outbam, "wb", template=samFile_list[0])
    if barcodes_out is None:
        barcodes_out = barcodes_in.copy()
    
    for ss in range(len(samFile_list)):
        samFile = samFile_list[ss]
        _barcodes_in = barcodes_in[ss]
        _barcodes_out = barcodes_out[ss]
        
        READ_CNT = 0
        reads_all = []
        for i in range(len(positions)):
            chrom = chroms[i]
            POS = positions[i]
            
            for _read in samFile.fetch(chrom, POS-1, POS):
                if _read.has_tag(cell_tag) == False:
                    continue
                try:
                    idx = _barcodes_in.index(_read.get_tag(cell_tag))
                    _read.set_tag(cell_tag, _barcodes_out[idx])
                except ValueError:
                    continue
                reads_all.append(_read)
                
                READ_CNT += 1
                if READ_CNT % 100000 == 0:
                    print("BAM%d: %.2fM reads." %(ss+1, READ_CNT/1000000))
        
        # remove redundant reads (one read may be called multiple times)
        reads_all = set(reads_all)
        print(len(reads_all), READ_CNT)
        for _read in reads_all:
            outbam.write(_read)
            
        samFile.close()
    outbam.close()
    return None


def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--samFiles", "-s", dest="sam_files", default=None,
        help=("Input bam or sam files, comma separated."))
    parser.add_option("--barcodeFiles", "-b", dest="barcodes_files", 
        default=None, help=("Input barcode files, comma separated."))
    parser.add_option("--regionFile", "-r", dest="region_file", 
        default=None, help=("Input SNP list."))
    parser.add_option("--doubletRate", "-d", dest="doublet_rate", 
        type="float", default=None, help=("Doublet rate [default: n/100000]"))
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Directory for output files: pooled.bam and barcodes_pool.tsv."))
    parser.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")

    group = OptionGroup(parser, "Cell barcodes sampling")
    group.add_option("--nCELL", type="int", dest="n_cell", default=None, 
        help="The number of cells in each sample [default: %default]")
    group.add_option("--minorSAMPLE", type="float", dest="minor_sample", 
        default=1.0, help="Ratio size of minor sample [default: %default]")
    group.add_option("--randomSEED", type="int", dest="random_seed", 
        default=None, help="The random seed in numpy [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to VCF_convert!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    ## out directory
    if options.out_dir is None:
        print("Error: need outDir for output files.")
        sys.exit(1)
    elif os.path.dirname(options.out_dir) == "":
        out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    ## sam files
    if options.sam_files is None:
        print("Error: need samFile for sam file.")
        sys.exit(1)
    else:
        samFile_list = options.sam_files.split(",")
    
    ## cell barcodes
    if options.barcodes_files is None:
        print("Error: need files for cell barcodes.")
        sys.exit(1)
    else:
        barcodes_files = options.barcodes_files.split(",")
    if len(barcodes_files) != len(samFile_list):
        print("Error: barcodes files are not equal to sam files.")
        sys.exit(1)
    
    barcodes_in = []
    for _bar in barcodes_files:
        fid = open(_bar, 'r')
        all_lines = [x.rstrip() for x in fid.readlines()]
        fid.close()
        barcodes_in.append(all_lines)
    if options.n_cell is not None:
        barcodes_in = sample_barcodes(barcodes_in, options.n_cell, 
            options.minor_sample, options.random_seed)
    barcodes_out = pool_barcodes(barcodes_in, out_dir, options.doublet_rate, 
        seed=options.random_seed)
    
    ## VCF file
    vcf_dat = load_VCF(options.region_file, biallelic_only=False, 
                       load_sample=False)
    chroms = vcf_dat['FixedINFO']['CHROM']
    positions = [int(x) for x in vcf_dat['FixedINFO']['POS']]
    
    # fetch each position
    if (options.nproc == 1):
        BAM_FILE = out_dir + "/pooled.bam"
        fetch_reads(samFile_list, chroms, positions, 
            BAM_FILE, barcodes_in, barcodes_out)
    else:
        result = []
        pool = multiprocessing.Pool(processes=options.nproc)
        for ii in range(len(samFile_list)):
            BAM_FILE = out_dir + "/pooled_temp%d.bam" %(ii)
            print(ii, BAM_FILE)
            result.append(pool.apply_async(fetch_reads, ([samFile_list[ii]], 
                chroms, positions, BAM_FILE, [barcodes_in[ii]], 
                [barcodes_out[ii]], "CB"), callback=show_progress))
        pool.close()
        pool.join()

        ## merge bam files
        file_list = [out_dir + "/pooled.bam"]
        file_list += [out_dir + "/pooled_temp%d.bam" %(x) 
                        for x in range(len(samFile_list))]
        bashCommand = "samtools merge %s" %(" ".join(file_list))
        pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        pro.communicate()[0]
        for dd in range(len(samFile_list)):
            os.remove(out_dir + "/pooled_temp%d.bam" %(dd))
    print("")

    ## sort and index bam file
    bashCommand = "samtools sort %s -o %s" %(out_dir + "/pooled.bam", 
        out_dir + "/pooled.sorted.bam")
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]
    
    bashCommand = "samtools index %s" %(out_dir + "/pooled.sorted.bam")
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]

    os.remove(out_dir + "/pooled.bam")
    
    
if __name__ == "__main__":
    main()
    
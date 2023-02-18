
import subprocess
import numpy as np
from scipy.io import mmread
from itertools import permutations

from .vcf_utils import load_VCF, write_VCF, parse_donor_GPb
from .vcf_utils import read_sparse_GeneINFO, GenoINFO_maker, match_SNPs

def match_donor_VCF(cell_dat, donor_vcf):
    """Match variants between cell VCF and donor VCF information
    """
    mm_idx = match_SNPs(cell_dat['variants'], donor_vcf['variants'])
    idx1 = np.where(mm_idx != None)[0] #remove None
    # TODO: check when chr is not compatible! given warning.
    if len(idx1) == 0:
        print("[vireo] warning: no variants matched to donor VCF, " + 
                "please check chr format!")
    else:
        print("[vireo] %d out %d variants matched to donor VCF" 
                %(len(idx1), len(cell_dat['variants'])))
    idx2 = mm_idx[idx1].astype(int)

    cell_dat['AD'] = cell_dat['AD'][idx1, :]
    cell_dat['DP'] = cell_dat['DP'][idx1, :]
    cell_dat["variants"]  = [cell_dat["variants"][x] for x in idx1]
    for _key in cell_dat["FixedINFO"].keys():
        cell_dat["FixedINFO"][_key] = [
            cell_dat["FixedINFO"][_key][x] for x in idx1]

    donor_vcf["variants"]  = [donor_vcf["variants"][x] for x in idx2]
    for _key in donor_vcf["FixedINFO"].keys():
        donor_vcf["FixedINFO"][_key] = [
            donor_vcf["FixedINFO"][_key][x] for x in idx2]
    for _key in donor_vcf["GenoINFO"].keys():
        donor_vcf["GenoINFO"][_key] = [
            donor_vcf["GenoINFO"][_key][x] for x in idx2]  

    return cell_dat, donor_vcf


def read_cellSNP(dir_name, layers=['AD', 'DP']):
    """Read data from the cellSNP output directory

    Parameters
    ----------
    dir_name:
        directory full path name for cellSNP output
    
    Return
    ------
    A disctionary containing AD, DP, cells and variants
    """
    cell_dat = load_VCF(dir_name + "/cellSNP.base.vcf.gz", load_sample=False,
                        biallelic_only=False)
    for _layer in layers:
        cell_dat[_layer] = mmread(dir_name + "/cellSNP.tag.%s.mtx" %(_layer)).tocsc()
    cell_dat['samples'] = np.genfromtxt(dir_name + "/cellSNP.samples.tsv", dtype=str)
    return cell_dat


def read_vartrix(alt_mtx, ref_mtx, cell_file, vcf_file=None):
    """Read data from VarTrix

    Parameters
    ----------
    alt_mtx:
        sparse matrix file for alternative alleles
    ref_mtx:
        sparse matrix file for reference alleles
    cell_file:
        file for cell barcodes, each per line
    vcf_file:
        the vcf file used for fetch variants in VarTrix
    
    Return
    ------
    A disctionary containing AD, DP, cells and optionally variants
    """
    if vcf_file is not None:
        cell_dat = load_VCF(vcf_file, load_sample=False, biallelic_only=False)
        cell_dat['variants'] = np.array(cell_dat['variants'])
    else:
        cell_dat = {}
    cell_dat['AD'] = mmread(alt_mtx).tocsc()
    cell_dat['DP'] = mmread(ref_mtx).tocsc() + cell_dat['AD']
    cell_dat['samples'] = np.genfromtxt(cell_file, dtype=str)
    return cell_dat


def write_donor_id(out_dir, donor_names, cell_names, n_vars, res_vireo):
    """
    Write the results of donor id into files.
    """
    ID_prob, doublet_prob = res_vireo['ID_prob'], res_vireo['doublet_prob']

    prob_max = np.max(ID_prob, axis=1)
    prob_doublet_out = np.max(doublet_prob, axis=1)
    donor_singlet = np.array(donor_names, "U100")[np.argmax(ID_prob, axis=1)]

    doublet_names = [",".join(x) for x in permutations(donor_names, 2)]
    donor_doublet = np.array(doublet_names, "U100")[np.argmax(doublet_prob, 
                                                              axis=1)]

    donor_ids = donor_singlet.copy()
    donor_ids[prob_max < 0.9] = "unassigned"
    donor_ids[prob_doublet_out >= 0.9] = "doublet"
    donor_ids[n_vars < 10] = "unassigned"


    ## save log file
    fid = open(out_dir + "/_log.txt", "w")
    fid.writelines("logLik: %.3e\n" %(res_vireo['LB_doublet']))
    fid.writelines("thetas: \n%s\n" %(res_vireo['theta_shapes']))
    fid.close()

    ## save summary file
    fid = open(out_dir + "/summary.tsv", "w")
    fid.writelines("Var1\tFreq\n")
    donor_ids_uniq, donor_ids_count = np.unique(donor_ids, return_counts=True)
    for i in range(len(donor_ids_uniq)):
        fid.writelines("%s\t%d\n" %(donor_ids_uniq[i], donor_ids_count[i]))
    fid.close()
    print("[vireo] final donor size:")
    print("\t".join([str(x) for x in donor_ids_uniq]))
    print("\t".join([str(x) for x in donor_ids_count]))

    ## save donor_ids file
    fid = open(out_dir + "/donor_ids.tsv", "w")
    header = ["cell", "donor_id", "prob_max", "prob_doublet", "n_vars", 
              "best_singlet", "best_doublet", "doublet_logLikRatio"]
    fid.writelines("\t".join(header) + "\n")
    for i in range(len(cell_names)):
        line = [cell_names[i], donor_ids[i], "%.2e" %prob_max[i],
                "%.2e" %prob_doublet_out[i], "%d" %n_vars[i],
                donor_singlet[i], donor_doublet[i], 
                "%.3f" %res_vireo['doublet_LLR'][i]]
        fid.writelines("\t".join(line) + "\n") 
    fid.close()

    ## save singlet probability file
    fid = open(out_dir + "/prob_singlet.tsv", "w")
    fid.writelines("\t".join(["cell"] + donor_names) + "\n")
    for i in range(len(cell_names)):
        line = ["%.2e" %x for x in ID_prob[i, :]]
        fid.writelines("\t".join([cell_names[i]] + line) + "\n") 
    fid.close()

    ## save doublet probability file
    fid = open(out_dir + "/prob_doublet.tsv", "w")
    fid.writelines("\t".join(["cell"] + doublet_names) + "\n")
    for i in range(len(cell_names)):
        line = ["%.2e" %x for x in doublet_prob[i, :]]
        fid.writelines("\t".join([cell_names[i]] + line) + "\n")
    fid.close()

    ## save ambient RNA file
    if res_vireo['ambient_Psi'] is not None:
        fid = open(out_dir + "/prop_ambient.tsv", "w")
        fid.writelines("\t".join(["cell"] + donor_names + ['logLik_ratio']) + "\n")
        for i in range(len(cell_names)):
            line = ["%.4e" %x for x in res_vireo['ambient_Psi'][i, :]]
            line += ['%.2f' %res_vireo['Psi_LLRatio'][i]]
            fid.writelines("\t".join([cell_names[i]] + line) + "\n")
        fid.close()

    bashCommand = "gzip -f %s %s" %(out_dir + "/prob_singlet.tsv", 
        out_dir + "/prob_doublet.tsv")
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]

def make_whitelists(donor_id_file, out_prefix):
    """Generate whitelist for each donor as input for umi_tools extract
    """
    table = np.genfromtxt(donor_id_file, dtype='str', delimiter='\t')[1:, :]
    table = table[table[:, 1] != 'unassigned', :]
    table = table[table[:, 1] != 'doublet', :]
    
    for _donor in np.unique(table[:, 1]):
        idx = table[:, 1] == _donor
        barcodes = table[idx, 0]
        fid = open(out_prefix + "_%s.txt" %_donor, "w")
        for _line in barcodes:
            fid.writelines(_line.split('-')[0] + '\n')
        fid.close()

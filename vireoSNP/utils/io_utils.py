
import subprocess
import numpy as np
from itertools import permutations

def write_donor_id(out_dir, donor_names, cell_names, n_vars, res_vireo):
    """
    Write the results of donor id into files.
    """
    ID_prob, doublet_prob = res_vireo['ID_prob'], res_vireo['doublet_prob']

    prob_max = np.max(ID_prob, axis=1)
    prob_doublet_out = np.sum(doublet_prob, axis=1)
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
    print("\t".join([str(x) for x in donor_ids_uniq]))
    print("\t".join([str(x) for x in donor_ids_count]))

    ## save donor_ids file
    fid = open(out_dir + "/donor_ids.tsv", "w")
    header = ["cell", "donor_id", "prob_max", "prob_doublet", "n_vars", 
              "best_singlet", "best_doublet"]
    fid.writelines("\t".join(header) + "\n")
    for i in range(len(cell_names)):
        line = [cell_names[i], donor_ids[i], "%.2e" %prob_max[i],
                "%.2e" %prob_doublet_out[i], "%d" %n_vars[i],
                donor_singlet[i], donor_doublet[i]]
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

    bashCommand = "gzip -f %s %s" %(out_dir + "/prob_singlet.tsv", 
        out_dir + "/prob_doublet.tsv")
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]

# def ppca_plot(AD, DP):
#     """
#     PPCA plot for each cell genotypes. This function is still underdevelopment
#     """
#     Z = DP.copy().astype(float)
#     idx = DP > 0
#     Z[idx] = AD[idx] / Z[idx]
#     Z[idx] = Z[idx] - 0.5

#     from sklearn.decomposition import TruncatedSVD
#     svd = TruncatedSVD(n_components=5, n_iter=7, random_state=42)
#     svd.fit(Z)

#     print("variance explained:", svd.explained_variance_ratio_)

#     import matplotlib.pyplot as plt
#     plt.scatter(svd.components_[0, :], svd.components_[1, :])
#     return svd.components_


def heat_matrix(X, yticks=None, xticks=None, rotation=45, cmap='BuGn', 
    alpha=0.6, **kwargs):
    """
    Plot heatmap of distance matrix
    """
    import matplotlib.pyplot as plt

    im = plt.imshow(X, cmap=cmap, alpha=alpha, **kwargs)
    plt.xticks(range(len(xticks)), xticks, rotation=rotation)
    plt.yticks(range(len(yticks)), yticks)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=rotation, ha="right",
    #          rotation_mode="anchor")
    
    # Loop over data dimensions and create text annotations.
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            plt.text(j, i, "%.2f" %X[i, j],
                     ha="center", va="center", color="k")
    
    return im

def plot_GT(out_dir, cell_GPb, donor_names, 
            donor_GPb=None, donor_names_in=None):
    """
    Plot the genotype distance between samples
    """
    import matplotlib.pyplot as plt

    ## compare the GT probability of estimated samples
    diff_mat = np.zeros((cell_GPb.shape[2], cell_GPb.shape[2]))
    for i in range(cell_GPb.shape[2]):
        for j in range(cell_GPb.shape[2]):
            diff_mat[i,j] = np.mean(np.abs(cell_GPb[:, :, i] - 
                                           cell_GPb[:, :, j]))

    fig = plt.figure()
    heat_matrix(diff_mat, donor_names, donor_names)
    plt.title("Geno Prob Delta: %d SNPs" %(cell_GPb.shape[0]))
    plt.tight_layout()
    fig.savefig(out_dir + "/fig_GT_distance_estimated.pdf", dpi=300)

    ## compare in the estimated sample with input samples
    if donor_GPb is not None:
        diff_mat = np.zeros((cell_GPb.shape[2], donor_GPb.shape[2]))
        for i in range(cell_GPb.shape[2]):
            for j in range(donor_GPb.shape[2]):
                diff_mat[i,j] = np.mean(np.abs(cell_GPb[:, :, i] - 
                                               donor_GPb[:, :, j]))

        fig = plt.figure()
        heat_matrix(diff_mat, donor_names, donor_names_in)
        plt.title("Geno Prob Delta: %d SNPs" %(cell_GPb.shape[0]))
        plt.tight_layout()
        fig.savefig(out_dir + "/fig_GT_distance_input.pdf", dpi=300)


def minicode_plot(barcode_set, var_ids=None, sample_ids=None, 
                  cmap="Set3"):
    import matplotlib.pyplot as plt
    
    mat = np.zeros((len(barcode_set[0][1:]), len(barcode_set)))
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            mat[i, j] = float(barcode_set[j][i + 1])
            
    im = plt.imshow(mat, cmap=cmap)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            plt.text(j, i, int(mat[i, j]), 
                     ha="center", va="center", color="k")
            
    if var_ids is None:
        var_ids = range(mat.shape[0])
    plt.yticks(range(mat.shape[0]), var_ids)
    
    if sample_ids is None:
        sample_ids = ["%s\nS%d" %(barcode_set[x], x)
                      for x in range(mat.shape[1])]
    else:
        sample_ids = ["%s\n%s" %(barcode_set[x], sample_ids[x])
                      for x in range(mat.shape[1])]
    plt.xticks(range(mat.shape[1]), sample_ids)
    
    return im
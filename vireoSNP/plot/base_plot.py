# base functions for plotting

import numpy as np


def heat_matrix(X, yticks=None, xticks=None, rotation=45, cmap='BuGn', 
                alpha=0.6, display_value=True, row_sort=False, 
                aspect='auto', **kwargs):
    """
    Plot heatmap of distance matrix
    """
    import matplotlib.pyplot as plt
    
    if row_sort:
        row_idx = np.argsort(np.dot(X, 2**np.arange(X.shape[1])))
        X = X[row_idx, :]

    im = plt.imshow(X, cmap=cmap, alpha=alpha, aspect=aspect, **kwargs)
    if xticks is not None:
        plt.xticks(range(len(xticks)), xticks, rotation=rotation)
    if yticks is not None:
        plt.yticks(range(len(yticks)), yticks)
    
    # Loop over data dimensions and create text annotations.
    if display_value:
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
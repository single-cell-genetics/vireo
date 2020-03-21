# base functions for plotting

import numpy as np

WeiZhu_colors = np.array(['#4796d7', '#f79e54', '#79a702', '#df5858', '#556cab', 
                          '#de7a1f', '#ffda5c', '#4b595c', '#6ab186', '#bddbcf', 
                          '#daad58', '#488a99', '#f79b78', '#ffba00'])

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
        plt.xlim(-0.5, len(xticks) - 0.5)
    if yticks is not None:
        plt.yticks(range(len(yticks)), yticks)
        plt.ylim(-0.5, len(yticks) - 0.5)
    
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
    diff_mat = np.zeros((cell_GPb.shape[1], cell_GPb.shape[1]))
    for i in range(cell_GPb.shape[1]):
        for j in range(cell_GPb.shape[1]):
            diff_mat[i,j] = np.mean(np.abs(cell_GPb[:, i, :] - 
                                           cell_GPb[:, j, :]))

    fig = plt.figure()
    heat_matrix(diff_mat, donor_names, donor_names)
    plt.title("Geno Prob Delta: %d SNPs" %(cell_GPb.shape[0]))
    plt.tight_layout()
    fig.savefig(out_dir + "/fig_GT_distance_estimated.pdf", dpi=300)

    ## compare in the estimated sample with input samples
    if donor_GPb is not None:
        diff_mat = np.zeros((cell_GPb.shape[1], donor_GPb.shape[1]))
        for i in range(cell_GPb.shape[1]):
            for j in range(donor_GPb.shape[1]):
                diff_mat[i,j] = np.mean(np.abs( cell_GPb[:, i, :] - 
                                               donor_GPb[:, j, :]))

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
    plt.yticks(range(len(var_ids)), var_ids)
    plt.ylim(-0.5, len(var_ids) - 0.5)
    
    if sample_ids is None:
        sample_ids = ["%s\nS%d" %(barcode_set[x], x)
                      for x in range(mat.shape[1])]
    else:
        sample_ids = ["%s\n%s" %(barcode_set[x], sample_ids[x])
                      for x in range(mat.shape[1])]
    plt.xticks(range(len(sample_ids)), sample_ids)
    plt.xlim(-0.5, len(sample_ids) - 0.5)

    return im

def anno_heat(X, anno, **kwargs):
    import seaborn as sns
    idx = np.argsort(np.dot(X, 2**np.arange(X.shape[1])) + 
                     anno * 2**X.shape[1])
    g = sns.clustermap(X[idx], cmap="GnBu", yticklabels=False,
                       col_cluster=False, row_cluster=False,
                       row_colors=WeiZhu_colors[anno][idx], **kwargs)
    
    for label in np.unique(anno):
        g.ax_col_dendrogram.bar(0, 0, color=WeiZhu_colors[label],
                                label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=6, title="True clone")
    g.cax.set_position([.95, .2, .03, .45])
    return g

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
# base functions for plotting

import numpy as np

vireo_colors = np.array(['#4796d7', '#f79e54', '#79a702', '#df5858', '#556cab', 
                         '#de7a1f', '#ffda5c', '#4b595c', '#6ab186', '#bddbcf', 
                         '#daad58', '#488a99', '#f79b78', '#ffba00'])

def heat_matrix(X, yticks=None, xticks=None, rotation=45, cmap='BuGn', 
                alpha=0.6, display_value=True, row_sort=False, 
                aspect='auto', interpolation='none', **kwargs):
    """
    Plot heatmap of distance matrix

    Parameters
    ----------
    X: numpy.array or matrix
        The matrix to plot in heatmap
    yticks: list
        The ticks ids for y axis
    xticks: list
        The ticks ids for x axis
    ratation: scalar
        The ratation angel for xticks
    cmap: str
        The colormap for the heatmap, more options: 
        https://matplotlib.org/stable/tutorials/colors/colormaps.html
    alpha: scalar
        The transparency, value between 0 and 1
    display_value: bool
        If True, dispaly the values in the heatmap
    raw_sort: bool
        If True, sort the rows with row index as
        row_idx = np.argsort(np.dot(X, 2**np.arange(X.shape[1])))
    aspect: str
        `aspect` in `plt.imshow`
    interpolation: str
        `interpolation` in `plt.imshow`
    **kwargs: keywords & values
        `**kwargs` for `plt.imshow`
    
    Returns
    -------
    The return from `plt.imshow`

    Examples
    --------

    .. plot::

        >>> from vireoSNP.plot import heat_matrix
        >>> import numpy as np
        >>> np.random.seed(1)
        >>> X = np.random.rand(5, 7)
        >>> heat_matrix(X)
    """
    import matplotlib.pyplot as plt
    
    if row_sort:
        row_idx = np.argsort(np.dot(X, 2**np.arange(X.shape[1])))
        X = X[row_idx, :]

    im = plt.imshow(X, cmap=cmap, alpha=alpha, aspect=aspect, 
                    interpolation=interpolation, **kwargs)
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
                  cmap="Set3", interpolation='none', **kwargs):
    import matplotlib.pyplot as plt
    
    mat = np.zeros((len(barcode_set[0][1:]), len(barcode_set)))
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            mat[i, j] = float(barcode_set[j][i + 1])
            
    im = plt.imshow(mat, cmap=cmap, interpolation=interpolation, **kwargs)
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


def anno_heat(X, row_anno=None, col_anno=None,
              row_order_ids=None, col_order_ids=None, 
              xticklabels=False, yticklabels=False,
              row_cluster=False, col_cluster=False,
              **kwargs):
    """
    Heatmap with column or row annotations. Based on seaborn.clustermap()
    Row or column will be ordered by the annotation group.
    
    Note, haven't tested if input both row_anno and col_anno.
    """
    
    import seaborn as sns
    
    # prepare row annotation
    if row_anno is not None:
        if row_order_ids is None:
            row_order_ids = list(np.unique(row_anno))
        else:
            row_order_ids = [x for x in row_order_ids]
        row_num = np.array([row_order_ids.index(x) for x in row_anno])

        dot_row = np.array(np.nansum(X, axis=1)).reshape(-1)
        idx_row = np.argsort(row_num * 2**X.shape[1])# + dot_row / dot_row.max())

        row_colors = vireo_colors[row_num][idx_row]
    else:
        row_colors = None
        row_order_ids = []
        idx_row = range(X.shape[0])
        
    # prepare col annotation
    if col_anno is not None:
        if col_order_ids is None:
            col_order_ids = list(np.unique(col_order_ids))
        else:
            col_order_ids = [x for x in col_order_ids]
        col_num = np.array([col_order_ids.index(x) for x in col_anno])

        dot_col = np.array(np.nansum(X, axis=0)).reshape(-1)
        idx_col = np.argsort(col_num * 2**X.shape[0])# + dot_row / dot_row.max())
        
        col_colors = vireo_colors[col_num][idx_col]
    else:
        col_colors = None
        col_order_ids = []
        idx_col = range(X.shape[1])
        
    ## plot with seaborn clustermap
    g = sns.clustermap(X[idx_row, :][:, idx_col], 
                       row_colors=row_colors, col_colors=col_colors,
                       col_cluster=col_cluster, row_cluster=row_cluster,
                       xticklabels=xticklabels, yticklabels=yticklabels,
                       **kwargs)
    
    if row_anno is not None:
        for i in range(len(row_order_ids)):
            g.ax_row_dendrogram.bar(0, 0, color=vireo_colors[i],
                                    label=row_order_ids[i], linewidth=0)
        g.ax_row_dendrogram.legend(loc="center", ncol=1, title="")
        
    if col_anno is not None:
        for i in range(len(col_order_ids)):
            g.ax_col_dendrogram.bar(0, 0, color=vireo_colors[i],
                                    label=col_order_ids[i], linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=6, title="")
    
    g.cax.set_position([1.01, .2, .03, .45])
    
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
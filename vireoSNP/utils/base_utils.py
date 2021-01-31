import numpy as np

def get_confusion(ids1, ids2):
    """Get confusion matrix
    
    Parameters
    ----------
    ids1: numpy.array or list
        id list in the first annotation
    ids2: numpy.array or list
        id list in the second annotation
        
    Return
    ------
    (confuse_mat, ids1_uniq, ids2_uniq)
    confuse_mat[i, j]: 
        number of samples have ids1 == ids1_uniq[i]
        and ids2 == id2_uniq[j]
    """
    if type(ids1) == list: ids1 = np.array(ids1)
    if type(ids2) == list: ids2 = np.array(ids2)
    
    ids1_uniq = np.unique(ids1)
    ids2_uniq = np.unique(ids2)
    
    confuse_mat = np.zeros((len(ids1_uniq), len(ids2_uniq)), dtype=int)
    for i, _id1 in enumerate(ids1_uniq):
        for j, _id2 in enumerate(ids2_uniq):
            confuse_mat[i, j] = np.sum((ids1 == _id1) * (ids2 == _id2))
            
    return confuse_mat, ids1_uniq, ids2_uniq

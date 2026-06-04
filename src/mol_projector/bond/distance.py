"""Pure numpy distance matrix computation."""

import numpy as np


def calc_distance_matrix(struc):
    """Calculate pairwise Euclidean distance matrix (non-periodic).

    Args:
        struc: Struc object

    Returns:
        (completed_distmtx, dist_mtx): For non-periodic systems, both are identical
        (shape [natom, natom]).
    """
    coords = struc.coord
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=-1))
    return dists, dists

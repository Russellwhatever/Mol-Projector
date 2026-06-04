"""Valence tracking for bond-order determination."""

import numpy as np
from .constants import CORRECT_VALENCE


def update_remain_valence(struc, bond_matrix=None):
    """Update valence, adj_degree, remain_valence for all atoms.

    Args:
        struc: Struc object
        bond_matrix: Optional bond matrix; if None, uses struc.bond_matrix
    """
    if bond_matrix is None:
        bond_matrix = struc.bond_matrix
    struc.adj_mtx = (bond_matrix != 0).astype(int)
    for i, atom in enumerate(struc.atom):
        atom.valence = np.sum(bond_matrix[i, :]).item()
        atom.adj_degree = np.sum(struc.adj_mtx[i, :]).item()
        if atom.ele in CORRECT_VALENCE:
            atom.remain_valence = CORRECT_VALENCE[atom.ele] - atom.valence
        else:
            atom.remain_valence = 9

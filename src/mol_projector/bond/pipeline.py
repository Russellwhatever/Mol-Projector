"""Main bond-order determination pipeline."""

import os
import numpy as np
from rdkit import Chem

from .valence import update_remain_valence
from .connectivity import generate_spanning_tree, find_rings, gen_multiple_bond, gen_multiple_bond_once
from .electron import modify_electron
from .molkit import struc2rdmol, save_fig


def get_bond_order(struc, output_smiles=False, visualize_path=None):
    """Determine bond orders for a molecular structure.

    This is the main orchestrator that applies the full bond-order determination pipeline:
    1. Try simple multi-bond generation
    2. Generate spanning tree (break cycles)
    3. Add bonds for large rings (5-7 membered)
    4. Add bonds for small rings (3-4 membered)
    5. Try multi-bond generation again
    6. Apply electron modifications (formal charge, carbene, conjugated systems)

    Args:
        struc: Struc object with 3D coordinates
        output_smiles: If True, also return SMILES string
        visualize_path: If set, save debug visualization at each step

    Returns:
        (status, rdmol, smiles_or_None):
            status: True (normal), '_e' (electron), '_cb' (carbene),
                    '_cs' (conjugated), False (failed)
            rdmol: RDKit Mol object with bond orders
            smiles_or_None: SMILES string if output_smiles is True
    """
    connectivity = struc.adj_mtx

    # Step 0: Check initial bond matrix
    if visualize_path:
        print("0. check original bond_matrix")
    update_remain_valence(struc)
    if visualize_path:
        from ..viz.structure_viz import visualize_structure
        visualize_structure(struc, visualize_path)

    # Step 1: Try simple multi-bond generation
    struc = gen_multiple_bond(struc, connectivity)
    update_remain_valence(struc)
    error = sum([atom.remain_valence if atom.ele in [1, 6, 7, 8, 15, 17, 35, 46] else 0
                 for atom in struc.atom])
    if error == 0:
        mol = struc2rdmol(struc)
        smiles = Chem.MolToSmiles(mol) if output_smiles else None
        return True, mol, smiles

    # Step 2: Spanning tree (remove cycles)
    if visualize_path:
        print("1. spanning tree")
    struc.bond_matrix = generate_spanning_tree(struc)
    update_remain_valence(struc)
    if visualize_path:
        from ..viz.structure_viz import visualize_structure
        visualize_structure(struc, visualize_path)

    # Step 3: Large rings (5-7 membered)
    if visualize_path:
        print("2. large rings")
    rings = find_rings(struc, connectivity, 5) + \
            find_rings(struc, connectivity, 6) + \
            find_rings(struc, connectivity, 7)
    for ring in rings:
        for index, i in enumerate(ring):
            i_ending_atom = ring[index - 1]
            if struc.atom[i].remain_valence and struc.atom[i_ending_atom].remain_valence:
                if not struc.bond_matrix[i][i_ending_atom]:
                    struc.bond_matrix[i][i_ending_atom] += 1
                    struc.bond_matrix[i_ending_atom][i] += 1
                    update_remain_valence(struc)
                    break
            elif struc.atom[i].remain_valence and struc.atom[index + 1].remain_valence:
                if not struc.bond_matrix[i][i_ending_atom]:
                    struc.bond_matrix[i][i_ending_atom] += 1
                    struc.bond_matrix[i_ending_atom][i] += 1
                    update_remain_valence(struc)
                    break
    if visualize_path:
        from ..viz.structure_viz import visualize_structure
        visualize_structure(struc, visualize_path)

    # Step 4: Small rings (3-4 membered)
    if visualize_path:
        print("4. smaller rings")
    rings = find_rings(struc, connectivity, 3) + \
            find_rings(struc, connectivity, 4)
    for ring in rings:
        for index, i in enumerate(ring):
            i_ending_atom = ring[index - 1]
            i_next_atom = ring[index + 1] if index < len(ring) - 1 else ring[0]
            if struc.atom[i].remain_valence and struc.atom[i_ending_atom].remain_valence:
                if not struc.bond_matrix[i][i_ending_atom]:
                    struc.bond_matrix[i][i_ending_atom] += 1
                    struc.bond_matrix[i_ending_atom][i] += 1
                    update_remain_valence(struc)
                    break
            elif struc.atom[i].remain_valence and struc.atom[i_next_atom].remain_valence:
                if not struc.bond_matrix[i][i_ending_atom]:
                    struc.bond_matrix[i][i_ending_atom] += 1
                    struc.bond_matrix[i_ending_atom][i] += 1
                    update_remain_valence(struc)
                    break
    if visualize_path:
        from ..viz.structure_viz import visualize_structure
        visualize_structure(struc, visualize_path)

    # Step 5: Multi-bond generation (second pass)
    if visualize_path:
        print("5. multiple bond")
    struc = gen_multiple_bond_once(struc, connectivity)
    update_remain_valence(struc)
    if visualize_path:
        from ..viz.structure_viz import visualize_structure
        visualize_structure(struc, visualize_path)

    # Step 6: Electron modifications
    if visualize_path:
        print("6 modify_electron")
    modify_result = modify_electron(struc, connectivity)
    if modify_result:
        if modify_result == "with electron":
            e_type = '_e'
        elif modify_result == 'carbene':
            e_type = '_cb'
        elif modify_result[0] == 'converge system':
            e_type = '_cs'
    else:
        e_type = ''
    if visualize_path:
        from ..viz.structure_viz import visualize_structure
        visualize_structure(struc, visualize_path)

    # Final check
    error = sum([atom.remain_valence if atom.ele in [1, 6, 7, 8, 15, 17, 35, 46] else 0
                 for atom in struc.atom])
    if error > 0:
        return False, None, None

    mol = struc2rdmol(struc)
    smiles = None
    if output_smiles:
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            pass
    return e_type if e_type else True, mol, smiles

"""Electron/charge modifications: formal charges, carbenes, conjugated systems."""

import numpy as np
from itertools import combinations
from .valence import update_remain_valence
from .constants import CORRECT_VALENCE


def modify_electron(struc, bond_matrix=None):
    """Handle remaining valence errors via formal charge, carbene, or conjugated systems.

    Returns:
        False: no modification made
        "with electron": formal charge adjustment applied
        "carbene": carbene (diradical) detected
        ("converge system", paths): conjugated pi-system bonding applied
    """
    if bond_matrix is None:
        bond_matrix = struc.bond_matrix

    def get_neighbors(s, i_atom, bond_order):
        neighbors = []
        for j in range(len(s.bond_matrix)):
            if bond_matrix[i_atom][j] > 0 and s.bond_matrix[i_atom][j] == bond_order:
                neighbors.append(j)
        return neighbors

    # Case 1: atoms with remain_valence == 2 (potential carbene or electron pair)
    valence_2_atoms = [i for i, atom in enumerate(struc.atom) if atom.remain_valence == 2]
    if valence_2_atoms:
        for i in valence_2_atoms:
            neighbors = []
            for j in [0, 1, 2]:
                neighbors.extend(get_neighbors(struc, i, j))
            carbene_maybe = True
            for neighbor in neighbors:
                if struc.atom[neighbor].remain_valence == 0:
                    if struc.atom[neighbor].ele in [7, 8]:
                        carbene_maybe = False
                        struc.bond_matrix[i][neighbor] += 1
                        struc.bond_matrix[neighbor][i] += 1
                        struc.atom[i].charge -= 1
                        struc.atom[neighbor].charge += 1
                        update_remain_valence(struc)
                        return "with electron"
            if carbene_maybe:
                struc.atom[i].remain_valence -= 2
                struc.atom[i].is_carbene = True
                return "carbene"

    # Case 2: Conjugated pi-systems (radical pairs with remain_valence == 1)
    def get_adjacent_atom_index(s, idx):
        return np.where(s.adj_mtx[idx, :] == 1)[0]

    def find_conjugate_pair(s):
        conjugate_pair_all = []
        for i, atom_i in enumerate(s.atom):
            if CORRECT_VALENCE[atom_i.ele] - atom_i.adj_degree > 0:
                adj_atoms = get_adjacent_atom_index(s, i)
                for j in adj_atoms:
                    if CORRECT_VALENCE[s.atom[j].ele] - s.atom[j].adj_degree > 0:
                        conjugate_pair = sorted([i, j])
                        if conjugate_pair not in conjugate_pair_all:
                            conjugate_pair_all.append(conjugate_pair)
        return conjugate_pair_all

    def find_all_paths(adj_mtx, start, end):
        if start == end:
            return [[start]]

        all_paths = []
        visited = set()

        def dfs(current, path):
            if current == end:
                all_paths.append(path[:])
                return

            visited.add(current)
            neighbors = np.where(adj_mtx[current] > 0)[0]

            for neighbor in neighbors:
                if neighbor not in visited:
                    dfs(neighbor, path + [neighbor])

            visited.remove(current)

        dfs(start, [start])
        return all_paths

    valence_1_atoms = [i for i, atom in enumerate(struc.atom)
                       if atom.remain_valence == 1 and CORRECT_VALENCE[atom.ele] - atom.adj_degree <= 1]
    if valence_1_atoms and len(valence_1_atoms) % 2 == 0:
        def generate_radical_pairs(atoms_list):
            if len(atoms_list) == 0:
                return []
            if len(atoms_list) == 2:
                return [[(atoms_list[0], atoms_list[1])]]

            pairs = []
            for i in range(1, len(atoms_list)):
                remaining = atoms_list[1:i] + atoms_list[i + 1:]
                sub_pairs = generate_radical_pairs(remaining)
                for sub_pair in sub_pairs:
                    pairs.append([(atoms_list[0], atoms_list[i])] + sub_pair)

            return pairs

        def is_valid_path(path):
            if len(path) < 3:
                return False

            middle_atoms = path[1:-1]

            for i, atom_idx in enumerate(middle_atoms):
                atom = struc.atom[atom_idx]

                if atom.remain_valence == 0 and atom.ele in [7, 8]:
                    continue

                if i > 0:
                    prev_atom_idx = middle_atoms[i - 1]
                    if sorted([atom_idx, prev_atom_idx]) in find_conjugate_pair(struc):
                        continue

                if i < len(middle_atoms) - 1:
                    next_atom_idx = middle_atoms[i + 1]
                    if sorted([atom_idx, next_atom_idx]) in find_conjugate_pair(struc):
                        continue

                return False

            return True

        all_pairings = generate_radical_pairs(valence_1_atoms)

        for pairing in all_pairings:
            all_paths_valid = True
            all_paths_found = []

            for atom1, atom2 in pairing:
                paths = find_all_paths(struc.adj_mtx, atom1, atom2)

                if not paths:
                    all_paths_valid = False
                    break

                valid_path_found = False
                for path in paths:
                    if is_valid_path(path):
                        all_paths_found.append(path)
                        valid_path_found = True
                        break

                if not valid_path_found:
                    all_paths_valid = False
                    break

            if all_paths_valid:
                for path in all_paths_found:
                    if len(path) % 2 == 0:
                        for i in range(0, len(path), 2):
                            if i + 1 < len(path):
                                struc.bond_matrix[path[i]][path[i + 1]] = 2
                                struc.bond_matrix[path[i + 1]][path[i]] = 2
                        for i in range(1, len(path) - 1, 2):
                            if i + 1 < len(path):
                                struc.bond_matrix[path[i]][path[i + 1]] = 1
                                struc.bond_matrix[path[i + 1]][path[i]] = 1
                    else:
                        for i in range(0, len(path) - 1, 2):
                            if i + 1 < len(path):
                                struc.bond_matrix[path[i]][path[i + 1]] = 2
                                struc.bond_matrix[path[i + 1]][path[i]] = 2
                        for i in range(1, len(path) - 1, 2):
                            if i + 1 < len(path):
                                struc.bond_matrix[path[i]][path[i + 1]] = 1
                                struc.bond_matrix[path[i + 1]][path[i]] = 1

                    update_remain_valence(struc)

                for path_atoms in all_paths_found:
                    if not sum(struc.atom[atom_idx].remain_valence for atom_idx in path_atoms):
                        for atom_idx in path_atoms:
                            struc.atom[atom_idx].charge -= struc.atom[atom_idx].remain_valence

                return "converge system", all_paths_found

    return False

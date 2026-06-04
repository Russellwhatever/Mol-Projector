"""Connectivity algorithms: spanning tree, ring finding, multi-bond generation."""

import numpy as np
from typing import List, Optional
from .valence import update_remain_valence
from .distance import calc_distance_matrix


def generate_spanning_tree(struc) -> np.ndarray:
    """Generate a minimum spanning tree bond matrix (no cycles).

    Uses Kruskal's algorithm with edge weights = interatomic distances.
    """
    def find(parent, i):
        if parent[i] != i:
            parent[i] = find(parent, parent[i])
        return parent[i]

    def union(parent, rank, x, y):
        root_x = find(parent, x)
        root_y = find(parent, y)
        if root_x == root_y:
            return False
        if rank[root_x] < rank[root_y]:
            parent[root_x] = root_y
        elif rank[root_x] > rank[root_y]:
            parent[root_y] = root_x
        else:
            parent[root_y] = root_x
            rank[root_x] += 1
        return True

    adj_matrix = (struc.bond_matrix != 0).astype(int)
    n_atoms = struc.natom

    edges = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if adj_matrix[i, j] > 0:
                dist = np.linalg.norm(struc.atom[i].xyz - struc.atom[j].xyz)
                edges.append((dist, i, j))

    edges.sort()
    parent = list(range(n_atoms))
    rank = [0] * n_atoms

    new_matrix = np.zeros((n_atoms, n_atoms))
    edge_count = 0
    for dist, i, j in edges:
        if union(parent, rank, i, j):
            new_matrix[i, j] = 1
            new_matrix[j, i] = 1
            edge_count += 1
            if edge_count == n_atoms - 1:
                break

    return new_matrix


def find_rings(struc, bond_matrix=None, ring_size: int = 5) -> list:
    """Find rings of a given size in the molecular structure using DFS.

    Args:
        struc: Struc object
        bond_matrix: Connectivity matrix; if None, uses struc.bond_matrix
        ring_size: Target ring size (number of atoms)

    Returns:
        List of rings, each ring is a list of atom indices
    """
    if bond_matrix is None:
        bond_matrix = struc.bond_matrix
    adj_matrix = (bond_matrix != 0).astype(int)
    n_atoms = struc.natom
    _, dist_matrix = calc_distance_matrix(struc)

    def is_valid_ring(path, check_every_distance=False):
        if not adj_matrix[path[-2]][path[0]]:
            return False
        if check_every_distance:
            for i in range(len(path)):
                for j in range(i + 1, len(path)):
                    if dist_matrix[path[i]][path[j]] > 2.6 * len(path) / 2:
                        return False
        return True

    def normalize(ring):
        unique = []
        for atom in ring:
            if atom not in unique:
                unique.append(atom)
        min_index = unique.index(min(ring))
        out_path = unique[min_index:] + unique[:min_index]
        return out_path

    def dfs(current, start, path, visited, rings):
        if len(path) >= ring_size + 1:
            if current == start:
                normalized_path = normalize(path)
                if normalized_path not in rings:
                    rings.append(normalized_path)
            return

        neighbors = np.where(adj_matrix[current] > 0)[0]
        for neighbor in neighbors:
            neighbor = int(neighbor)
            if (neighbor != start and not visited[neighbor]) or \
               (neighbor == start and len(path) >= ring_size):
                visited[neighbor] = True
                path.append(neighbor)
                dfs(neighbor, start, path, visited, rings)
                path.pop()
                visited[neighbor] = False

    all_rings = []
    for start_atom in range(n_atoms):
        visited = [False] * n_atoms
        visited[start_atom] = True
        path = [start_atom]
        dfs(start_atom, start_atom, path, visited, all_rings)

    rings = []
    sorted_rings = []
    for ring in all_rings:
        if sorted(ring) not in sorted_rings:
            sorted_rings.append(sorted(ring))
            rings.append(ring)

    rings.sort(key=lambda x: (len(x), x))
    return rings


def gen_multiple_bond(in_struc, bond_matrix: Optional[np.ndarray] = None):
    """Generate multiple bonds iteratively: N/O first (2 passes), then C (2 passes).

    Modifies struc in-place.
    """
    if bond_matrix is None:
        bond_matrix = in_struc.bond_matrix

    def get_adjacent_atom_index(bm, idx):
        adjacent_indices = np.where(bm[idx, :] == 1)[0]
        return sorted(adjacent_indices, key=lambda i: in_struc.gen_dist(idx, i))

    def gen(struc, O_or_N: bool):
        for i, atom_i in enumerate(struc.atom):
            if O_or_N:
                if atom_i.ele > 6 and atom_i.remain_valence:
                    adjacent_indices = get_adjacent_atom_index(bond_matrix, i)
                    for j in adjacent_indices:
                        if struc.atom[j].remain_valence and struc.bond_matrix[i][j] > 0:
                            struc.bond_matrix[i][j] += 1
                            struc.bond_matrix[j][i] += 1
                            update_remain_valence(struc)
                            break
            else:
                if atom_i.ele == 6 and atom_i.remain_valence:
                    adjacent_indices = get_adjacent_atom_index(bond_matrix, i)
                    for j in adjacent_indices:
                        if struc.atom[j].remain_valence and struc.bond_matrix[i][j] > 0:
                            struc.bond_matrix[i][j] += 1
                            struc.bond_matrix[j][i] += 1
                            update_remain_valence(struc)
                            break

    gen(in_struc, O_or_N=True)
    gen(in_struc, O_or_N=True)
    gen(in_struc, O_or_N=False)
    gen(in_struc, O_or_N=False)

    return in_struc


def gen_multiple_bond_once(in_struc, bond_matrix: Optional[np.ndarray] = None):
    """Generate multiple bonds in a single pass (all elements together, 2 passes).

    Modifies struc in-place.
    """
    if bond_matrix is None:
        bond_matrix = in_struc.bond_matrix

    def get_adjacent_atom_index(bm, idx):
        adjacent_indices = np.where(bm[idx, :] == 1)[0]
        return sorted(adjacent_indices, key=lambda i: in_struc.gen_dist(idx, i))

    def gen(struc):
        for i, atom_i in enumerate(struc.atom):
            if atom_i.remain_valence:
                adjacent_indices = get_adjacent_atom_index(bond_matrix, i)
                for j in adjacent_indices:
                    if struc.atom[j].remain_valence and struc.bond_matrix[i][j] > 0:
                        struc.bond_matrix[i][j] += 1
                        struc.bond_matrix[j][i] += 1
                        update_remain_valence(struc)
                        break

    gen(in_struc)
    gen(in_struc)
    return in_struc

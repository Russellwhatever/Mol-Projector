"""Tests for bond-order determination."""

import numpy as np
from mol_projector.core import Struc
from mol_projector.bond import update_remain_valence, calc_distance_matrix
from mol_projector.bond.connectivity import generate_spanning_tree, gen_multiple_bond_once
from mol_projector.bond.pipeline import get_bond_order


def make_water():
    s = Struc()
    s.add_atom(8, [0, 0, 0.117])
    s.add_atom(1, [0, 0.757, -0.469])
    s.add_atom(1, [0, -0.757, -0.469])
    return s


def make_formaldehyde():
    s = Struc()
    s.add_atom(6, [-0.526, 0.0, 0.0])
    s.add_atom(8, [-0.526, -1.215, 0.0])
    s.add_atom(1, [0.574, 0.0, 0.0])
    s.add_atom(1, [-1.099, 0.938, 0.0])
    return s


class TestValence:
    def test_water(self):
        s = make_water()
        update_remain_valence(s)
        # O: valence=2, remain=0; H: valence=1, remain=0
        assert s.atom[0].remain_valence == 0  # O satisfied
        assert s.atom[1].remain_valence == 0  # H satisfied
        assert s.atom[2].remain_valence == 0  # H satisfied

    def test_update_after_bond_change(self):
        s = make_water()
        update_remain_valence(s)
        s.bond_matrix[0, 1] += 1  # upgrade one O-H to double
        s.bond_matrix[1, 0] += 1
        update_remain_valence(s)
        assert s.atom[0].remain_valence < 0  # O over-bonded


class TestDistance:
    def test_simple(self):
        s = Struc()
        s.add_atom(6, [0, 0, 0])
        s.add_atom(8, [3, 4, 0])
        completed, dist = calc_distance_matrix(s)
        assert abs(dist[0, 1] - 5.0) < 0.001
        assert abs(dist[1, 0] - 5.0) < 0.001
        assert dist[0, 0] == 0


class TestSpanningTree:
    def test_water(self):
        s = make_water()
        tree = generate_spanning_tree(s)
        # Tree should have n-1 = 2 edges, no cycles
        assert np.sum(tree) // 2 == 2
        # Check acyclic: sum of vertex degrees in tree
        degrees = np.sum(tree, axis=1)
        assert np.all(degrees > 0)


class TestPipeline:
    def test_water(self):
        s = make_water()
        status, mol, smiles = get_bond_order(s, output_smiles=True)
        assert status == True
        assert smiles is not None
        assert 'O' in smiles

    def test_formaldehyde(self):
        s = make_formaldehyde()
        status, mol, smiles = get_bond_order(s, output_smiles=True)
        assert status == True
        # Should have C=O double bond
        assert s.bond_matrix[0, 1] == 2  # C=O

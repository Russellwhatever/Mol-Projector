"""Tests for core data model."""

import numpy as np
from mol_projector.core import Atom, Struc, copy_struc
from mol_projector.core.periodictable import ele_dict, ele_table, ele_mass


class TestPeriodicTable:
    def test_ele_dict(self):
        assert ele_dict['H'] == 1
        assert ele_dict['C'] == 6
        assert ele_dict['O'] == 8
        assert ele_dict['Pd'] == 46

    def test_ele_table(self):
        assert ele_table[0] == 'H'
        assert ele_table[5] == 'C'
        assert ele_table[7] == 'O'

    def test_ele_mass(self):
        assert len(ele_mass) >= 80
        assert ele_mass[0] == 1.008


class TestAtom:
    def test_create(self):
        a = Atom(6, [1.0, 2.0, 3.0])
        assert a.ele == 6
        assert np.allclose(a.xyz, [1.0, 2.0, 3.0])
        assert a.charge == 0

    def test_ele_symb(self):
        a = Atom(6, [0, 0, 0])
        assert a.ele_symb == 'C'

    def test_set_ele_symb(self):
        a = Atom(1, [0, 0, 0])
        a.ele_symb = 'O'
        assert a.ele == 8
        assert a.ele_symb == 'O'

    def test_valence_attrs(self):
        a = Atom(6, [0, 0, 0])
        assert a.valence == 0
        assert a.adj_degree == 0
        assert a.remain_valence == 0
        assert not a.is_carbene

    def test_hash(self):
        a1 = Atom(6, [1.0, 2.0, 3.0])
        a2 = Atom(6, [1.0, 2.0, 3.0])
        assert hash(a1) != hash(a2)  # different objects


class TestStruc:
    def test_create_empty(self):
        s = Struc()
        assert s.natom == 0

    def test_add_atom(self):
        s = Struc()
        s.add_atom(8, [0, 0, 0.117])
        s.add_atom(1, [0, 0.757, -0.469])
        s.add_atom(1, [0, -0.757, -0.469])
        assert s.natom == 3
        assert list(s.iza) == [8, 1, 1]

    def test_computed_properties(self):
        s = Struc()
        s.add_atom(6, [0, 0, 0])
        s.add_atom(8, [1.2, 0, 0])
        assert s.natom == 2
        assert s.coord.shape == (2, 3)
        assert list(s.iza) == [6, 8]
        assert s.ele_list == [6, 8]

    def test_bond_matrix_auto_init(self):
        s = Struc()
        s.add_atom(8, [0, 0, 0])
        s.add_atom(1, [0, 0.96, 0])
        s.add_atom(1, [0, -0.96, 0])
        bm = s.bond_matrix
        assert bm.shape == (3, 3)
        assert bm[0, 1] >= 1  # O-H bond
        assert bm[0, 2] >= 1  # O-H bond

    def test_gen_dist(self):
        s = Struc()
        s.add_atom(6, [0, 0, 0])
        s.add_atom(8, [1.2, 0, 0])
        assert abs(s.gen_dist(0, 1) - 1.2) < 0.001

    def test_copy_struc(self):
        s = Struc()
        s.add_atom(6, [0, 0, 0])
        s.add_atom(8, [1.2, 0, 0])
        s.molecule_name = 'test'
        _ = s.bond_matrix  # trigger init
        s2 = copy_struc(s)
        assert s2.natom == 2
        assert s2.molecule_name == 'test'
        assert s2.bond_matrix.shape == (2, 2)

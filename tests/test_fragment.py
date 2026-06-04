"""Tests for molecular fragment decomposition."""

import os
import tempfile
from mol_projector.io import read_structure
from mol_projector.bond import get_bond_order
from mol_projector.fragment import report_fragments


def test_water_no_surrounding():
    struc = read_structure(
        os.path.join(os.path.dirname(__file__), 'data', 'water.arc'), fmt='arc'
    )[0]
    get_bond_order(struc)

    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
        path = f.name

    try:
        bond_info = report_fragments(
            struc, save_path=path, bicyclic=True, double_bond=True,
            contain_surrounding_atom=False
        )
        # O-H bonds should be found
        assert len(bond_info) == 2
        assert [0, 1] in bond_info  # O-H
        assert [0, 2] in bond_info  # O-H

        with open(path) as f:
            content = f.read()
        assert content.startswith('1\n')  # 1 fragment
    finally:
        os.unlink(path)


def test_formaldehyde_no_surrounding():
    struc = read_structure(
        os.path.join(os.path.dirname(__file__), 'data', 'formaldehyde.arc'), fmt='arc'
    )[0]
    get_bond_order(struc)

    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
        path = f.name

    try:
        bond_info = report_fragments(
            struc, save_path=path, bicyclic=True, double_bond=True,
            contain_surrounding_atom=False
        )
        # C=O double bond merges the two fragments into one
        assert len(bond_info) == 3  # C=O, C-H, C-H
    finally:
        os.unlink(path)


def test_formaldehyde_with_surrounding():
    struc = read_structure(
        os.path.join(os.path.dirname(__file__), 'data', 'formaldehyde.arc'), fmt='arc'
    )[0]
    get_bond_order(struc)

    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
        path = f.name

    try:
        bond_info = report_fragments(
            struc, save_path=path, bicyclic=True, double_bond=True,
            contain_surrounding_atom=True
        )
        # With surrounding atoms, should still find all bonds
        assert len(bond_info) >= 3
    finally:
        os.unlink(path)

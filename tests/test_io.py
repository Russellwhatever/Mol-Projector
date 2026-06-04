"""Tests for file I/O."""

import os
from mol_projector.io import read_structure, list_supported_formats, ReadPath


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class TestReadPath:
    def test_scan_arc(self):
        rp = ReadPath(DATA_DIR, fmt='arc')
        assert 'water' in rp
        assert 'formaldehyde' in rp

    def test_file2name(self):
        rp = ReadPath(DATA_DIR, fmt='arc')
        name = rp.file2name('/some/path/molecule.arc')
        assert name == 'molecule'


class TestARCReader:
    def test_read_water(self):
        path = os.path.join(DATA_DIR, 'water.arc')
        struc_list = read_structure(path, fmt='arc')
        assert len(struc_list) == 1
        s = struc_list[0]
        assert s.natom == 3
        assert list(s.iza) == [8, 1, 1]

    def test_read_formaldehyde(self):
        path = os.path.join(DATA_DIR, 'formaldehyde.arc')
        struc_list = read_structure(path, fmt='arc')
        assert len(struc_list) == 1
        s = struc_list[0]
        assert s.natom == 4
        assert 6 in s.ele_list  # C
        assert 8 in s.ele_list  # O


class TestFormats:
    def test_supported_formats(self):
        fmts = list_supported_formats()
        assert 'arc' in fmts
        assert 'xyz' in fmts
        assert 'mol' in fmts
        assert 'sdf' in fmts

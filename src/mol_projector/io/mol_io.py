"""MOL/SDF format readers and writers via RDKit."""

import os
import numpy as np
from rdkit import Chem
from . import register_reader
from ..bond.molkit import rdmol2struc, save_fig


@register_reader('mol')
def read_mol(filename: str):
    """Read a MOL file via RDKit.

    Args:
        filename: Path to the .mol file

    Returns:
        List of Struc objects (typically one)
    """
    mol = Chem.MolFromMolFile(filename, removeHs=False)
    if mol is None:
        print(f"Warning: could not read {filename}")
        return []
    struc = rdmol2struc(mol)
    struc.molecule_name = os.path.splitext(os.path.basename(filename))[0]
    return [struc]


@register_reader('sdf')
def read_sdf(filename: str):
    """Read an SDF file via RDKit.

    Args:
        filename: Path to the .sdf file

    Returns:
        List of Struc objects
    """
    supplier = Chem.SDMolSupplier(filename, removeHs=False)
    result = []
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        struc = rdmol2struc(mol)
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"{i}"
        struc.molecule_name = name
        result.append(struc)
    return result


def write_mol(mol, path: str):
    """Write an RDKit Mol to a MOL file.

    Args:
        mol: RDKit Mol object
        path: Output file path
    """
    Chem.MolToMolFile(mol, path)


def write_mol_with_fig(mol, mol_path: str, fig_path: str = None):
    """Write MOL file and optionally a 2D depiction figure.

    Args:
        mol: RDKit Mol object
        mol_path: Path for the .mol output file
        fig_path: Optional path for .png figure output
    """
    os.makedirs(os.path.dirname(mol_path), exist_ok=True)
    Chem.MolToMolFile(mol, mol_path)
    if fig_path:
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        save_fig(mol, fig_path)

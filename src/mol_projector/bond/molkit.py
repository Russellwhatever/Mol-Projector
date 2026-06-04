"""RDKit interoperability: Struc <-> Mol conversion, figure generation."""

from rdkit import Chem
from ..core.struc import copy_struc


def struc2rdmol(struc, remove_H=False, only_single_bond=False):
    """Convert a Struc object to an RDKit Mol object.

    Args:
        struc: Struc object with bond_matrix
        remove_H: If True, remove hydrogen atoms
        only_single_bond: If True, all bonds are single bonds

    Returns:
        RDKit Chem.Mol object with 3D conformer
    """
    struc_new = copy_struc(struc)

    if remove_H:
        new_atoms = [atom for atom in struc_new.atom if atom.ele != 1]
        struc_new.atom = new_atoms

    atomic_nums = [atom.ele for atom in struc_new.atom]
    coordinates = struc_new.coord.tolist()
    natom = struc_new.natom
    bond_matrix = struc_new.bond_matrix

    mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    for idx, atom in enumerate(struc_new.atom):
        rdkit_atom = Chem.Atom(atom.ele)
        if hasattr(atom, 'charge') and atom.charge != 0:
            rdkit_atom.SetFormalCharge(int(round(atom.charge)))
        if hasattr(atom, 'is_carbene') and atom.is_carbene:
            rdkit_atom.SetNumRadicalElectrons(2)
        mol.AddAtom(rdkit_atom)

    for i in range(natom):
        for j in range(i):
            bond_order = bond_matrix[i][j]
            if not only_single_bond:
                if atomic_nums[i] == 46 and atomic_nums[j] == 15 and bond_order:
                    mol.AddBond(j, i, Chem.rdchem.BondType.DATIVE)
                    continue
                if atomic_nums[j] == 46 and atomic_nums[i] == 15 and bond_order:
                    mol.AddBond(i, j, Chem.rdchem.BondType.DATIVE)
                    continue
                if bond_order == 1:
                    mol.AddBond(i, j, Chem.rdchem.BondType.SINGLE)
                elif bond_order == 2:
                    mol.AddBond(i, j, Chem.rdchem.BondType.DOUBLE)
                elif bond_order == 3:
                    mol.AddBond(i, j, Chem.rdchem.BondType.TRIPLE)
                elif bond_order == 4:
                    mol.AddBond(i, j, Chem.rdchem.BondType.SINGLE)
            else:
                if bond_order > 0:
                    mol.AddBond(i, j, Chem.rdchem.BondType.SINGLE)

    mol = mol.GetMol()
    conf = Chem.rdchem.Conformer(natom)
    for idx, xyz in enumerate(coordinates):
        conf.SetAtomPosition(idx, xyz)
    mol.AddConformer(conf)
    return mol


def rdmol2struc(rdmol) -> "Struc":
    """Convert an RDKit Mol object to a Struc object."""
    from ..core.struc import Struc

    def get_pos_from_mol(mol):
        if not mol.GetNumConformers():
            return [[0.0, 0.0, 0.0] for _ in range(mol.GetNumAtoms())]
        conformer = mol.GetConformer()
        atom_coords = []
        for i in range(mol.GetNumAtoms()):
            try:
                pos = conformer.GetAtomPosition(i)
                coord = [pos.x, pos.y, pos.z]
            except Exception:
                coord = [0.0, 0.0, 0.0]
            atom_coords.append(coord)
        return atom_coords

    def get_ele_from_mol(mol):
        return [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    def get_bond_matrix_from_mol(mol):
        import numpy as np
        num_atom = mol.GetNumAtoms()
        bond_matrix = np.zeros((num_atom, num_atom), dtype=int)
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            order = int(bond_type) if bond_type != Chem.rdchem.BondType.DATIVE else 1
            bond_matrix[i, j] = order
            bond_matrix[j, i] = order
        return bond_matrix

    ele = get_ele_from_mol(rdmol)
    pos = get_pos_from_mol(rdmol)
    struc = Struc()
    struc.max_f = 0.0

    for i in range(len(ele)):
        struc.add_atom(ele[i], pos[i])
    struc.bond_matrix = get_bond_matrix_from_mol(rdmol)

    import numpy as np
    max_coord = np.max(struc.coord, axis=0)
    min_coord = np.min(struc.coord, axis=0)
    pbc_span = max_coord - min_coord
    pbc_span_max = np.max(pbc_span) + 10
    struc.abc = [pbc_span_max, pbc_span_max, pbc_span_max, 90, 90, 90]
    return struc


def save_fig(mol, filename, remove_H=True):
    """Save a 2D depiction of the molecule as a PNG file."""
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit import Chem

    Chem.rdDepictor.Compute2DCoords(mol, clearConfs=True)
    if remove_H:
        mol = Chem.rdmolops.RemoveHs(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    drawer.WriteDrawingText(filename)

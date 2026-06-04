"""XYZ format parser."""

import numpy as np
from . import register_reader
from ..core import Struc
from ..core.periodictable import ele_dict


@register_reader('xyz')
def read_xyz(filename: str):
    """Parse an XYZ file.

    Standard XYZ format:
        Line 1: number of atoms
        Line 2: comment line (may contain energy or lattice info)
        Lines 3+: element_symbol x y z

    Args:
        filename: Path to the .xyz file

    Returns:
        List of Struc objects (typically one)
    """
    with open(filename) as f:
        lines = f.readlines()

    if not lines:
        return []

    result = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        try:
            natom = int(line)
        except ValueError:
            i += 1
            continue

        comment = lines[i + 1].strip() if i + 1 < len(lines) else ""

        struc = Struc()
        # Try to extract energy from comment line
        if 'energy=' in comment.lower():
            try:
                energy_str = comment.split('energy=')[1].split()[0]
                struc.energy = float(energy_str)
            except Exception:
                pass

        # Parse lattice info if present (ASE-style extended XYZ)
        if 'Lattice=' in comment:
            try:
                lattice_str = comment.split('Lattice="')[1].split('"')[0]
                lattice = np.array([float(x) for x in lattice_str.split()])
                if len(lattice) == 9:
                    struc.abc = struc.lat2abc(lattice.reshape(3, 3))
            except Exception:
                pass

        for j in range(natom):
            if i + 2 + j >= len(lines):
                break
            parts = lines[i + 2 + j].strip().split()
            if len(parts) < 4:
                continue
            elem_sym = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            ele = ele_dict.get(elem_sym)
            if ele is None:
                continue
            struc.add_atom(ele, [x, y, z])

        if struc.natom > 0:
            struc.molecule_name = os.path.basename(filename)
            result.append(struc)

        i += 2 + natom

    return result

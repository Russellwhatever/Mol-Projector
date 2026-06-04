import math as m
import hashlib
import numpy as np
from typing import List, Union, Literal
from . import periodictable as pt
from .atom import Atom

list_or_array = Union[List[int], np.ndarray]

COVALENT_RADIUS = [
    0.25, 1.20,
    1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 1.60,
    1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.70,
    2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.80,
    2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.90,
    2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.80, 1.60, 1.90, 2.00,
]

MAX_BOND_NUMBER = {
    1: 1, 4: 6, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1,
    12: 6, 13: 6, 14: 6, 15: 6, 16: 4, 17: 1, 35: 1, 46: 2,
}

COMBINE2CUTOFF = {
    2: ((0.85, 1),),
    42: ((1.2, 1),),
    56: ((1.15, 1),),
    72: ((1.1, 1),),
    870: ((1.8, 4),),
    432: ((1.25, 3), (1.43, 2), (1.64, 1)),
    546: ((1.26, 3), (1.4, 2), (1.6, 1)),
    672: ((1.15, 3), (1.30, 2), (1.55, 1)),
    6090: ((2.3, 4),),
    686: ((1.2, 3), (1.35, 2), (1.55, 1)),
    840: ((1.15, 3), (1.3, 2), (1.56, 1)),
    7308: ((2.35, 4),),
    1024: ((1.3, 2), (1.58, 1)),
    8584: ((2.4, 4),),
    48778: ((2.8, 5),),
}


class Struc(object):
    def __init__(self):
        super(Struc, self).__init__()

        self.atom: List[Atom] = []
        self.abc: np.ndarray = np.zeros(6)
        self.energy: float = 0
        self._max_f: float = 0

        self.label: str = None
        self.molecule_name: str = None

        self._bond_matrix: np.ndarray = None
        self._adj_mtx: np.ndarray = None

    def __getattribute__(self, __name):
        if __name == "natom":
            return len(self.atom)
        elif __name == "lat":
            return self.abc2lat(self.abc)
        elif __name == "coord":
            return np.array([at.xyz for at in self.atom])
        elif __name == "iza":
            return np.array([at.ele for at in self.atom])
        elif __name == "charge":
            return np.array([at.charge for at in self.atom])
        elif __name == "force":
            return np.array([at.force for at in self.atom])
        elif __name == "max_f":
            if np.max(abs(self.force)) == 0:
                return self._max_f
            else:
                return np.max(abs(self.force))
        elif __name == "ele_compos":
            _ary = np.array(np.unique(self.iza, return_counts=True)).T
            return {k: v for k, v in _ary}
        elif __name == "ele_list":
            return list(self.ele_compos.keys())
        elif __name == "ele_natom":
            return list(self.ele_compos.values())
        elif __name == "ele":
            return [atom.ele for atom in self.atom]
        elif __name == "nele":
            return len(self.ele_compos)
        elif __name == "ele_name_list" or __name == "ele_symb_list":
            return [pt.ele_table[ele - 1] for ele in sorted(self.ele_compos)]
        elif __name == "centroid":
            return np.sum(self.coord, axis=0) / len(self.atom)
        elif __name == "coord1D":
            return self.coord.reshape(-1)
        elif __name == "fcoord":
            return np.matmul(self.coord, np.linalg.inv(self.lat))
        elif __name == "fcoord1D":
            return self.fcoord.reshape(-1)
        elif __name == "short_repr":
            if hasattr(self, "formula"):
                name = self.formula
            elif hasattr(self, "smiles_name"):
                name = self.smiles_name
            elif hasattr(self, "ecfp_name"):
                name = self.ecfp_name
            else:
                name = str(self.ele_compos)
            return name
        else:
            return super().__getattribute__(__name)

    def __getattr__(self, __name):
        if __name == 'bond_matrix':
            result = self._calc_bond_matrix_dist()
            result = self._cut_bond_matrix(result)
            self.__dict__['bond_matrix'] = result
            return result
        elif __name == 'adj_mtx':
            result = (self.bond_matrix != 0).astype(int)
            self.__dict__['adj_mtx'] = result
            return result
        else:
            return super().__getattr__(__name)

    def __setattr__(self, __name, __value):
        if __name == "lat":
            self.abc = self.lat2abc(__value)
        elif __name == "coord":
            for at, xyz in zip(self.atom, __value):
                at.xyz = xyz
        elif __name == "iza":
            for at, ele in zip(self.atom, __value):
                at.ele = ele
        elif __name == "charge":
            for at, charge in zip(self.atom, __value):
                at.charge = charge
        elif __name == "force":
            for at, force in zip(self.atom, __value):
                at.force = force
        elif __name == "fcoord":
            self.coord = np.matmul(__value, self.lat)
        elif __name == "max_f":
            self._max_f = __value
        elif __name == "bond_matrix":
            self._bond_matrix = __value
            self.__dict__['bond_matrix'] = __value
        elif __name == "adj_mtx":
            self._adj_mtx = __value
            self.__dict__['adj_mtx'] = __value
        else:
            return super(Struc, self).__setattr__(__name, __value)

    def __hash__(self):
        _tmp = hashlib.md5()
        for attr in [self.abc, self.energy, self.natom]:
            _tmp.update(str(attr).encode())
        rng = np.random.RandomState(89)
        indx = rng.randint(low=0, high=self.natom, size=1000)
        for attr in [self.iza, self.charge]:
            _tmp.update(attr.flat[indx].data)
        indx = rng.randint(low=0, high=self.natom * 3, size=1000)
        for attr in [self.coord, self.force]:
            _tmp.update(attr.flat[indx].data)
        return _tmp.hexdigest()

    def __repr__(self):
        return "<%s obj %.40s at %s>" % (__class__.__name__, self.short_repr,
                                         id(self))

    def abc2lat(self, in_abc: list_or_array) -> np.ndarray:
        a, b, c = in_abc[0:3]
        alpha, beta, gamma = [x * np.pi / 180.0 for x in in_abc[3:]]
        bc2 = b**2 + c**2 - 2 * b * c * m.cos(alpha)
        h1 = a
        h2 = b * m.cos(gamma)
        h3 = b * m.sin(gamma)
        h4 = c * m.cos(beta)
        h5 = ((h2 - h4)**2 + h3**2 + c**2 - h4**2 - bc2) / (2 * h3)
        h6 = m.sqrt(c**2 - h4**2 - h5**2)
        return np.array([[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]])

    def lat2abc(self, in_lat: np.ndarray) -> np.ndarray:
        lat = in_lat
        a = np.linalg.norm(lat[0])
        b = np.linalg.norm(lat[1])
        c = np.linalg.norm(lat[2])
        alpha = m.acos(np.dot(lat[1], lat[2]) / (b * c)) * 180.0 / np.pi
        beta = m.acos(np.dot(lat[0], lat[2]) / (a * c)) * 180.0 / np.pi
        gamma = m.acos(np.dot(lat[0], lat[1]) / (a * b)) * 180.0 / np.pi
        return np.array([a, b, c, alpha, beta, gamma])

    def add_atom(self,
                 ele: int,
                 xyz: list_or_array,
                 charge: float = .0,
                 force: list_or_array = None,
                 **kwargs) -> None:
        self.atom.append(Atom(ele=ele, xyz=xyz, charge=charge, force=force, **kwargs))

    def add_force(self, force: list_or_array, index: int) -> None:
        if type(force) == list or type(force) == np.ndarray:
            self.atom[index].force = force
        elif type(force) == str:
            self.atom[index].force = [float(x) for x in force.split()]

    def add_stress(self, stress: Union[np.ndarray, str]) -> None:
        if type(stress) == np.ndarray:
            self.stress = stress
        elif type(stress) == str:
            float_stress = [float(x) for x in stress.split()]
            assert float_stress.__len__() == 6, ValueError("wrong float information")
            self.stress = float_stress

    def sort_atom(self,
                  key: Literal["ele", "x", "y", "z"] = "ele") -> np.ndarray:
        if key == "ele":
            idx = np.argsort(self.iza)
        elif key == "x":
            idx = np.argsort(self.coord[:, 0])
        elif key == "y":
            idx = np.argsort(self.coord[:, 1])
        elif key == "z":
            idx = np.argsort(self.coord[:, 2])
        self.atom = [self.atom[i] for i in idx]
        return idx

    def get_direction(self, iatom1: int, iatom2: int) -> np.ndarray:
        f_vec = np.array(self.fcoord[iatom1]) - np.array(self.fcoord[iatom2])
        vec = np.matmul(np.array([x - np.round(x) for x in f_vec]), self.lat)
        return vec

    def gen_dist(self, iatom1: int, iatom2: int) -> float:
        return np.linalg.norm(self.atom[iatom1].xyz - self.atom[iatom2].xyz)

    def gen_dist_mtx(self) -> np.ndarray:
        n = self.natom
        dist_mtx = np.zeros([n, n])
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(self.atom[i].xyz - self.atom[j].xyz)
                dist_mtx[i, j] = d
                dist_mtx[j, i] = d
        return dist_mtx

    def gen_mass_centre(self):
        centre = np.array([0., 0., 0.])
        mass = 0
        for atom in self.atom:
            centre += atom.xyz / pt.ele_mass[atom.ele]
            mass += pt.ele_mass[atom.ele]
        self.mass_centre = centre / mass
        return self.mass_centre

    def translation(self, trans_xyz: np.ndarray = np.array([0, 0, 0])):
        if not hasattr(self, "mass_centre"):
            self.gen_mass_centre()
        for atom in self.atom:
            atom.xyz += trans_xyz
        self.mass_centre += trans_xyz

    def rotation(self, theta: float = 0, phi: float = 0):
        if not hasattr(self, "mass_centre"):
            self.gen_mass_centre()
        xyz = np.array([atom.xyz - self.mass_centre for atom in self.atom]).T
        rot_mtx_xy = np.array([[np.cos(theta), np.sin(theta), 0],
                               [-np.sin(theta), np.cos(theta), 0],
                               [0, 0, 1]])
        rot_mtx_yz = np.array([[1, 0, 0],
                               [0, np.cos(phi), np.sin(phi)],
                               [0, -np.sin(phi), np.cos(phi)]])
        rot_mtx = np.dot(rot_mtx_xy, rot_mtx_yz)
        xyz = np.dot(rot_mtx, xyz).T
        for i, atom in enumerate(self.atom):
            atom.xyz = xyz[i] + self.mass_centre
        return rot_mtx

    def delete_atom(self, ele: int, xyz: np.ndarray):
        for iatom, atom in enumerate(self.atom):
            if atom.ele == ele and not (atom.xyz - xyz).any():
                self.atom.pop(iatom)
                return True
        return False

    def expand_cell(self, a: int, b: int, c: int):
        old_atom = self.atom
        self.atom = []
        for ia in range(a):
            for ib in range(b):
                for ic in range(c):
                    delta_coord = np.matmul(self.lat.T, np.array([ia, ib, ic]))
                    for atom in old_atom:
                        self.add_atom(atom.ele, atom.xyz + delta_coord)
        self.abc[:3] *= np.array([a, b, c])

    def _calc_bond_matrix_dist(self, scale_factor: float = 1.0) -> np.ndarray:
        """Calculate initial bond matrix using distance cutoffs (pure numpy)."""
        n = self.natom
        bm = np.zeros((n, n), dtype=int)
        coords = self.coord
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(coords[i] - coords[j])
                ele_i, ele_j = self.atom[i].ele, self.atom[j].ele
                key = ele_i * ele_j * (ele_i + ele_j)
                if key in COMBINE2CUTOFF:
                    for cutoff, max_order in COMBINE2CUTOFF[key]:
                        if d <= cutoff * scale_factor:
                            bm[i, j] = max_order
                            bm[j, i] = max_order
                            break
                else:
                    ri = COVALENT_RADIUS[ele_i] if ele_i < len(COVALENT_RADIUS) else 1.5
                    rj = COVALENT_RADIUS[ele_j] if ele_j < len(COVALENT_RADIUS) else 1.5
                    if d <= (ri + rj) * 1.2 * scale_factor:
                        bm[i, j] = 1
                        bm[j, i] = 1
        return bm

    def _cut_bond_matrix(self, input_matrix: np.ndarray) -> np.ndarray:
        """Limit the number of bonds per atom using MAX_BOND_NUMBER."""
        output_matrix = (input_matrix != 0).astype(int)
        dist_mtx = self.gen_dist_mtx()

        for i, atom_i in enumerate(self.atom):
            max_bond = MAX_BOND_NUMBER.get(atom_i.ele, 9)
            if np.sum(output_matrix[i, :]) > max_bond:
                bonded_atoms = np.where(output_matrix[i, :] > 0)[0]
                bond_orders = [(j, input_matrix[i, j]) for j in bonded_atoms]
                bond_distances = [(j, dist_mtx[i, j]) for j in bonded_atoms]
                sorted_pairs = sorted(zip(bond_distances, bond_orders), key=lambda x: x[0][1])
                bond_distances, bond_orders = zip(*sorted_pairs)
                expanded_bond_orders = []
                for m, n in bond_orders:
                    for _ in range(n):
                        expanded_bond_orders.append((m, n))
                for j, _ in expanded_bond_orders[max_bond:]:
                    output_matrix[i, j] = 0
                    output_matrix[j, i] = 0

        return output_matrix

    def calc_bond_info(self) -> List[dict]:
        bond_matrix = self.bond_matrix
        adjacent_matrix = bond_matrix != 0
        ele_all = np.array([atom.ele for atom in self.atom])
        ele_list = sorted(self.ele_list)
        ele_mask = {ele: (ele_all == ele).repeat(self.natom).reshape(-1, self.natom).transpose() for ele in ele_list}
        connect_num = {ele: np.sum(mask * adjacent_matrix, axis=1) for ele, mask in ele_mask.items()}
        connect_num_list = [{ele: connect_num[ele][idx] for ele in ele_list} for idx in range(self.natom)]
        return connect_num_list


def copy_struc(struc: Struc) -> Struc:
    new_struc = Struc()
    new_struc.molecule_name = getattr(struc, 'molecule_name')
    new_struc.label = getattr(struc, 'molecule_name')
    new_struc.energy = getattr(struc, 'energy')
    new_struc.max_f = getattr(struc, 'max_f')
    new_struc.abc = getattr(struc, 'abc')

    if 'bond_matrix' in struc.__dict__:
        new_struc.__dict__['bond_matrix'] = np.copy(struc.__dict__['bond_matrix'])

    for atom in struc.atom:
        new_struc.add_atom(atom.ele, atom.xyz, atom.charge)
    return new_struc

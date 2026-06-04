import numpy as np
from typing import Union, List

from . import periodictable as pt

list_or_array = Union[List[int], np.ndarray]


class Atom(object):
    def __init__(self,
                 ele: int,
                 xyz: list_or_array,
                 charge: int = 0,
                 force: list_or_array = None,
                 tag: int = -1):
        self.ele = ele
        self.xyz = np.array(xyz, dtype=float)
        self.charge = charge
        self.force = np.array(force if force is not None else [0, 0, 0], dtype=float)
        self.is_carbene = False
        self.tag = tag

        self.valence: float = 0.0
        self.adj_degree: int = 0
        self.remain_valence: int = 0

    def __setattr__(self, __name, __value):
        if __name == "ele_symb":
            self.ele = pt.ele_dict[__value]
        else:
            return super().__setattr__(__name, __value)

    def __getattribute__(self, __name):
        if __name == "ele_symb":
            return pt.ele_table[self.ele - 1]
        else:
            return super().__getattribute__(__name)

    def __getattr__(self, __name):
        if __name == 'fix':
            self.fix = [False, False, False]
            return self.fix
        else:
            return super().__getattr__(__name)

    def __repr__(self):
        return "<Atom %s_%.3f_%.3f_%.3f in struc_%s>" % (
            self.ele_symb, self.xyz[0], self.xyz[1], self.xyz[2],
            id(self))

    def __hash__(self) -> int:
        return hash((id(self), float(self.xyz[0]), float(self.xyz[1]),
                     float(self.xyz[2]), self.ele))

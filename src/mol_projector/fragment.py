"""分子片段拆分模块。

将分子按环结构拆分为片段，支持桥环/螺环合并、双键连接合并、包含周边原子。
"""

import os
import numpy as np
from typing import List, Tuple

from mol_projector.core.struc import Struc
from mol_projector.bond.connectivity import find_rings


def report_fragments(
    struc: Struc,
    save_path: str = "./report",
    bicyclic: bool = True,
    double_bond: bool = True,
    contain_surrounding_atom: bool = True,
) -> List[Tuple[int, int]]:
    """将分子按环结构拆分为片段并输出报告文件。

    算法：
    1. 查找 3/4/5/6/7 元环，去重（子环被父环包含则移除）
    2. bicyclic=True: 合并共享原子的桥环/螺环
    3. 环上原子附加 H；环外非 H 重原子单独成片段
    4. double_bond=True: 双键/三键连接的片段合并
    5. contain_surrounding_atom=True: 每个片段额外包含一个邻接重原子

    Args:
        struc: 含 bond_matrix 的结构
        save_path: 报告输出路径
        bicyclic: 合并桥环/螺环
        double_bond: 合并多重键连接的片段
        contain_surrounding_atom: 包含相邻重原子

    Returns:
        片段内键连关系列表 [(atom_i, atom_j), ...]，索引从 0 开始
    """
    if struc is None or struc.natom == 0:
        print("Warning: Invalid structure input")
        return []

    save_dir = os.path.dirname(save_path)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    # ---- 查找环 ----
    rings = []
    for ring_range in [3, 4, 5, 6, 7]:
        rings.extend(find_rings(struc, None, ring_range))

    def _remove_duplicate_rings(rings_list: list) -> list:
        if not rings_list:
            return []
        # 先按排序后的内容去重（相同原子集合只保留一个）
        seen = set()
        deduped = []
        for r in rings_list:
            key = tuple(sorted(r))
            if key not in seen:
                seen.add(key)
                deduped.append(r)
        # 再移除真子集（小环被大环完全包含则移除）
        rings_list = sorted(deduped, key=len)
        keep = [True] * len(rings_list)
        for i, ring_i in enumerate(rings_list):
            set_i = set(ring_i)
            for j, ring_j in enumerate(rings_list):
                if i != j and set_i < set(ring_j):
                    keep[i] = False
                    break
        return [r for r, k in zip(rings_list, keep) if k]

    def _append_H(atom_idx: int) -> list:
        return [int(i) for i in np.where(struc.bond_matrix[atom_idx] > 0)[0]
                if struc.atom[int(i)].ele == 1]

    rings = _remove_duplicate_rings(rings)

    # ---- 1. 合并桥环/螺环 ----
    if bicyclic and rings:
        merged = True
        while merged:
            merged = False
            new_rings = []
            used = set()
            for i, current_ring in enumerate(rings):
                if i in used:
                    continue
                merged_ring = current_ring.copy()
                used.add(i)
                for j, other_ring in enumerate(rings):
                    if j in used:
                        continue
                    if set(current_ring) & set(other_ring):
                        for atom in other_ring:
                            if atom not in merged_ring:
                                merged_ring.append(atom)
                        used.add(j)
                        merged = True
                new_rings.append(merged_ring)
            rings = new_rings
        rings = _remove_duplicate_rings(rings)

    # ---- 构建片段 ----
    fragments = []
    atoms_in_rings = set()
    for ring in rings:
        atoms_in_rings.update(ring)
        for atom in list(ring):
            ring.extend(_append_H(atom))
        fragments.append(ring)

    for i in range(struc.natom):
        if i not in atoms_in_rings and struc.atom[i].ele != 1:
            frag = [i]
            frag.extend(_append_H(i))
            fragments.append(frag)

    # ---- 2. 多重键连接合并 ----
    def _get_ring_from_atom(atom_idx: int, in_frags: list) -> list:
        for frag in in_frags:
            if atom_idx in frag:
                return frag
        return []

    if double_bond and fragments:
        for fragment in fragments:
            if not fragment:
                continue
            checked = set()
            for atom_i in fragment:
                if atom_i in checked:
                    continue
                checked.add(atom_i)
                for conn in np.where(struc.bond_matrix[atom_i] > 1)[0]:
                    if conn not in fragment:
                        fragments.append(fragment + _get_ring_from_atom(conn, fragments))

    fragments = _remove_duplicate_rings(fragments)

    # ---- 3. 包含周边重原子 ----
    def _get_surrounding_atoms(ring: list) -> list:
        result = []
        for atom in ring:
            for nb in np.where(struc.bond_matrix[atom] > 0)[0]:
                nb = int(nb)
                if struc.atom[nb].ele != 1 and nb not in ring:
                    result.append(nb)
        return result

    frags_with_surroundings = []
    for fragment in fragments:
        frags_with_surroundings.append(fragment + _get_surrounding_atoms(fragment))

    fragments = _remove_duplicate_rings(fragments)
    frags_with_surroundings = _remove_duplicate_rings(frags_with_surroundings)

    # ---- 片段内键连关系 ----
    def _bond_pairs(frags_list: list) -> list:
        if not frags_list:
            return []
        pairs = []
        for frag in frags_list:
            if not frag or len(frag) < 2:
                continue
            for idx, a in enumerate(frag):
                for b in frag[idx + 1:]:
                    if (b, a) not in pairs and struc.bond_matrix[a][b] > 0:
                        pairs.append([a, b])
        return sorted(pairs)

    # ---- 写文件 ----
    def _write(frags_list: list, path: str):
        with open(path, 'w') as fo:
            fo.write(f"{len(frags_list)}\n")
            for frag in frags_list:
                if frag and len(frag) > 1:
                    fo.write(f"{len(frag)}\n")
                    fo.write(f"{' '.join(str(a + 1) for a in frag)}\n")
        print(f"Report saved to: {path}")

    if contain_surrounding_atom:
        _write(sorted(frags_with_surroundings, reverse=True, key=len), save_path)
        return _bond_pairs(frags_with_surroundings)
    else:
        _write(sorted(fragments, reverse=True, key=len), save_path)
        return _bond_pairs(fragments)


def fragment_cli(args=None):
    """片段拆分的命令行入口。"""
    import argparse

    parser = argparse.ArgumentParser(
        description='将分子按环结构拆分为片段',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input', type=str, help='输入文件路径')
    parser.add_argument('-f', '--format', type=str, default='arc',
                        choices=['arc', 'mol', 'xyz', 'sdf'], help='输入格式')
    parser.add_argument('-o', '--output', type=str, default='./report',
                        help='输出路径 (default: ./report)')
    parser.add_argument('--bicyclic', action='store_true', default=True,
                        help='合并桥环/螺环 (default)')
    parser.add_argument('--no-bicyclic', action='store_false', dest='bicyclic')
    parser.add_argument('--double-bond', action='store_true', default=True,
                        help='合并双键/三键连接片段 (default)')
    parser.add_argument('--no-double-bond', action='store_false', dest='double_bond')
    parser.add_argument('--surrounding', action='store_true', default=True,
                        help='包含相邻重原子 (default)')
    parser.add_argument('--no-surrounding', action='store_false', dest='surrounding')

    args = parser.parse_args(args)

    # 读取
    from mol_projector.io import read_structure
    from mol_projector.bond import get_bond_order
    from mol_projector.bond.molkit import rdmol2struc

    if args.format == 'mol':
        from rdkit import Chem
        mol = Chem.MolFromMolFile(args.input, removeHs=False)
        if mol is None:
            raise SystemExit(f"Failed to read: {args.input}")
        struc = rdmol2struc(mol)
    else:
        slist = read_structure(args.input, fmt=args.format)
        if not slist:
            raise SystemExit(f"No structures in: {args.input}")
        struc = slist[0]

    # 判定键级
    get_bond_order(struc)

    # 拆分
    bond_info = report_fragments(
        struc,
        save_path=args.output,
        bicyclic=args.bicyclic,
        double_bond=args.double_bond,
        contain_surrounding_atom=args.surrounding,
    )

    # 输出键连关系
    if bond_info:
        unique = sorted({(min(a, b), max(a, b)) for a, b in bond_info})
        bond_path = args.output + ".bond"
        with open(bond_path, 'w') as f:
            for a, b in unique:
                f.write(f"{a + 1} {b + 1}\n")
        print(f"Bond info saved to: {bond_path}")


if __name__ == '__main__':
    fragment_cli()

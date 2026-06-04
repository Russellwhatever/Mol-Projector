# Python API

## 核心数据类型

### Atom — 原子

```python
from mol_projector.core import Atom

atom = Atom(ele=6, xyz=[0.0, 0.0, 0.0], charge=0)

# 属性
atom.ele              # 6                 — 原子序数
atom.ele_symb         # 'C'               — 元素符号（自动查表）
atom.xyz              # np.array([0,0,0]) — 笛卡尔坐标
atom.charge           # 0                 — 形式电荷
atom.is_carbene       # False             — 是否为卡宾（双自由基）
atom.valence          # 0.0               — 当前键级总和
atom.adj_degree       # 0                 — 邻接原子数
atom.remain_valence   # 0                 — 剩余化合价
```

### Struc — 分子结构

```python
from mol_projector.core import Struc

struc = Struc()
struc.add_atom(8, [0, 0, 0.117])          # 添加 O 原子
struc.add_atom(1, [0, 0.757, -0.469])     # 添加 H 原子
struc.add_atom(1, [0, -0.757, -0.469])    # 添加 H 原子

# 计算属性（首次访问时自动计算）
struc.natom           # 3            — 原子总数
struc.coord           # (3,3) 数组   — 所有原子坐标
struc.iza             # [8, 1, 1]    — 所有原子序数
struc.ele_symb_list   # ['H', 'O']   — 元素种类（去重排序）
struc.ele_compos      # {8:1, 1:2}   — 元素组成
struc.lat             # (3,3) 数组   — 晶格向量
struc.fcoord          # (3,3) 数组   — 分数坐标

# bond_matrix — 首次访问时根据距离自动计算键连矩阵
struc.bond_matrix     # (3,3) 数组   — 键级矩阵
struc.adj_mtx         # (3,3) 数组   — 邻接矩阵（键级>0 → 1）
```

### 辅助函数

```python
from mol_projector.core import copy_struc

# 深拷贝结构（含键连矩阵）
s2 = copy_struc(struc)
```

## 文件读写

```python
from mol_projector.io import read_structure, ReadPath, list_supported_formats

# 查看所有支持的格式
list_supported_formats()  # ['arc', 'xyz', 'mol', 'sdf']

# 扫描目录 → {文件名: 文件路径}
rp = ReadPath('/data/molecules', fmt='arc')
for name, path in rp.items():
    print(f'{name} → {path}')

# 读取单个文件 → 返回 Struc 列表（多帧文件如 SDF 可能有多个）
struc_list = read_structure('water.arc', fmt='arc')
struc = struc_list[0]   # 取第一个结构
```

## 键级判定

```python
from mol_projector.bond import get_bond_order

status, rdmol, smiles = get_bond_order(struc, output_smiles=True)
```

**返回值说明：**

| `status` | 含义 |
|----------|------|
| `True` | 正常判定成功 |
| `'_e'` | 通过电子 / 电荷修正成功 |
| `'_cb'` | 通过卡宾（双自由基）修正成功 |
| `'_cs'` | 通过共轭体系修正成功 |
| `False` | 判定失败，存在无法消除的化合价错误 |

**rdmol** 是 RDKit `Chem.Mol` 对象，可自由使用 RDKit 的全部功能：

```python
from rdkit import Chem

# 写入 MOL 文件
Chem.MolToMolFile(rdmol, 'output.mol')

# 生成 SMILES
smiles = Chem.MolToSmiles(rdmol)

# 读取键级信息
for bond in rdmol.GetBonds():
    print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
```

## 输出 MOL 和结构图

```python
from mol_projector.io.mol_io import write_mol_with_fig

# 同时输出 MOL 和 2D 结构图
write_mol_with_fig(rdmol, 'output.mol', 'output.png')
```

## 键级相关工具函数

```python
from mol_projector.bond import update_remain_valence, calc_distance_matrix
from mol_projector.bond.molkit import struc2rdmol, rdmol2struc, save_fig

# 更新化合价状态
update_remain_valence(struc)

# 计算非周期性距离矩阵
completed, dist = calc_distance_matrix(struc)

# Struc ↔ RDKit Mol 互转
mol = struc2rdmol(struc)
struc2 = rdmol2struc(mol)

# 单独保存 2D 结构图
save_fig(mol, 'molecule.png')
```

## 可视化模块（可选）

```python
from mol_projector.viz import visualize_structure

# 在调试管线的每一步调用，保存当前键连状态图
visualize_structure(struc, 'step_output.png')
# 输出：节点=原子（序号+元素），边=键连（灰=单键、蓝=双键、红=三键）
```

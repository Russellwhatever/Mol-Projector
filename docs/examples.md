# 示例

## 水分子 H₂O

**输入（water.arc）**：

```
 Energy 0.0 0.0
 O      0.000000000    0.000000000    0.117000000  CORE  0  0  0  0  0
 H      0.000000000    0.757000000   -0.469000000  CORE  0  0  0  0  0
 H      0.000000000   -0.757000000   -0.469000000  CORE  0  0  0  0  0
 end
```

**Python 调用**：

```python
from mol_projector.io import read_structure
from mol_projector.bond import get_bond_order

struc = read_structure('water.arc', fmt='arc')[0]
status, mol, smiles = get_bond_order(struc, output_smiles=True)

print(f'status={status}, smiles={smiles}')
# status=True, smiles=[H]O[H]

print(struc.bond_matrix)
# [[0 1 1]
#  [1 0 0]
#  [1 0 0]]
```

**分析**：距离初始化后 O 与两个 H 各距离约 0.96 Å，在 O-H 截断距离（1.1 Å）以内，形成两个单键。O 化合价 = 2✓，H 化合价 = 1✓，直接成功。

---

## 甲醛 H₂CO

**输入（formaldehyde.arc）**：

```
 Energy 0.0 0.0
 C     -0.526000000    0.000000000    0.000000000  CORE  0  0  0  0  0
 O     -0.526000000   -1.215000000    0.000000000  CORE  0  0  0  0  0
 H      0.574000000    0.000000000    0.000000000  CORE  0  0  0  0  0
 H     -1.099000000    0.938000000    0.000000000  CORE  0  0  0  0  0
 end
```

**Python 调用**：

```python
struc = read_structure('formaldehyde.arc', fmt='arc')[0]
status, mol, smiles = get_bond_order(struc, output_smiles=True)

print(f'status={status}, smiles={smiles}')
# status=True, smiles=[H]C([H])=O

print(struc.bond_matrix)
# [[0 2 1 1]    ← C=O 双键, C-H 单键 × 2
#  [2 0 0 0]
#  [1 0 0 0]
#  [1 0 0 0]]
```

**管线追踪**：

| 步骤 | 操作 | 中间状态 |
|------|------|----------|
| 0 | 距离初始化 | C-O 距离 1.215 Å → 落入 C-O 键的双键范围（1.15~1.30 Å），但初始仅标记连通性 |
| 1 | 多重键尝试 | C 剩 1 价，O 剩 1 价 → `gen_multiple_bond` 将 C-O 升级为双键 |
| 检查 | 化合价验证 | C=4 ✓、O=2 ✓、H=1 ✓ → 直接成功 |

---

## 分子片段拆分

`mol-fragment` 命令将分子按环结构拆分为片段。

```bash
# 基本用法
mol-fragment molecule.arc -f arc -o ./frag_report

# 不合并双键连接（每个环/重原子独立成片段）
mol-fragment molecule.arc --no-double-bond

# 不包含周边重原子（输出纯片段）
mol-fragment molecule.arc --no-surrounding
```

**Python API**：

```python
from mol_projector.fragment import report_fragments
from mol_projector.io import read_structure
from mol_projector.bond import get_bond_order

struc = read_structure('molecule.arc', fmt='arc')[0]
get_bond_order(struc)
bond_info = report_fragments(struc, save_path='./report', bicyclic=True,
                              double_bond=True, contain_surrounding_atom=False)
```

输出格式：

```
2              ← 片段数
4              ← 片段 0 原子数
1 3 4 5 2      ← 片段 0 原子序号（1-based）
3              ← 片段 1 原子数
6 7 8          ← 片段 1 原子序号
```

同时生成 `.bond` 文件记录片段内键连关系。

---

## 批量转换

```python
from mol_projector.io import ReadPath, read_structure
from mol_projector.bond import get_bond_order
from mol_projector.io.mol_io import write_mol_with_fig
import os

name2path = ReadPath('/data/molecules', fmt='arc')
os.makedirs('output/mol', exist_ok=True)
os.makedirs('output/fig', exist_ok=True)

for name, path in name2path.items():
    struc_list = read_structure(path, fmt='arc')
    status, mol, smiles = get_bond_order(struc_list[0], output_smiles=True)

    if status:
        mol_path = f'output/mol/{name}.mol'
        fig_path = f'output/fig/{name}.png'
        write_mol_with_fig(mol, mol_path, fig_path)
        print(f'{name}: 成功 (status={status}, SMILES={smiles})')
    else:
        print(f'{name}: 失败')
```

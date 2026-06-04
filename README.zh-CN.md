# mol-projector

[![English](https://img.shields.io/badge/README-English-blue)](README.md)

将三维分子结构转换为二维 MOL 格式，并自动判定键级（单键/双键/三键）。

对于仅给出原子三维坐标的分子，`mol-projector` 利用化合价规则自动判定键级，处理特殊情况（卡宾、共轭体系、形式电荷），输出带有正确键级信息的 MOL 文件。

## 安装

```bash
# 从源码安装
cd mol-projector
pip install -e .

# 使用 conda
conda env create -f env.yaml
conda activate mol-projector
```

详见 [docs/install.md](docs/install.md)。

## 命令行用法

```bash
# 转换目录中所有 .arc 文件
mol-projector -i /path/to/arc/files -o /path/to/output

# 开启调试可视化
mol-projector -i . -o ./converted -d ./debug_vis

# 指定输入格式
mol-projector -i . -o ./converted --format xyz
```

详见 [docs/usage.md](docs/usage.md)（完整参数表、输出结构、报告解读）。

## Python API

```python
from mol_projector.io import read_structure
from mol_projector.bond import get_bond_order

# 读取结构文件
struc_list = read_structure('molecule.arc', fmt='arc')
struc = struc_list[0]

# 判定键级
status, rdmol, smiles = get_bond_order(struc, output_smiles=True)
# status 含义: True (正常), '_e' (电子修正), '_cb' (卡宾修正),
#            '_cs' (共轭体系修正), False (失败)

# 输出结果
from mol_projector.io.mol_io import write_mol_with_fig
write_mol_with_fig(rdmol, 'output.mol', 'output.png')
```

详见 [docs/api.md](docs/api.md)（完整模块接口、数据类型说明）。

## 算法流程

详见 [docs/algorithm.md](docs/algorithm.md)（七步管线详解、调试可视化）。

键级判定管线分为以下步骤：

1. **初始键连矩阵** — 根据原子间距离和元素对截断半径（`COMBINE2CUTOFF`）生成初始连接关系
2. **多重键生成** — 根据化合价需求将单键升级为双键/三键
3. **生成树** — 打破环结构，建立基础连接骨架
4. **大环闭合** — 重新添加 5-7 元环的键
5. **小环闭合** — 重新添加 3-4 元环的键
6. **多重键补全** — 第二轮双键/三键分配
7. **电子修正** — 处理形式电荷、卡宾（双自由基）以及共轭 π 体系

## 支持的格式

| 格式 | 扩展名 | 状态 |
|--------|-----------|--------|
| ARC (BIOVIA/LASP) | `.arc` | 已支持 |
| XYZ | `.xyz` | 已支持 |
| MOL (MDL Molfile) | `.mol` | 已支持 |
| SDF | `.sdf` | 已支持 |

详见 [docs/formats.md](docs/formats.md)（格式详解、示例、如何扩展新格式）。

## 测试

```bash
pip install pytest
pytest tests/ -v
```

## 验证案例

| 分子 | 输入格式 | SMILES | 键级 |
|----------|-------|--------|-------------|
| 水 (H2O) | ARC | `[H]O[H]` | O-H 单键 |
| 甲醛 (H2CO) | ARC | `[H]C([H])=O` | C=O 双键 |

详见 [docs/examples.md](docs/examples.md)（批量转换示例、管线过程追踪）。

---

## 文档目录

| 文档 | 说明 |
|------|------|
| [docs/install.md](docs/install.md) | 安装指南 — pip / conda 安装、环境配置、验证 |
| [docs/usage.md](docs/usage.md) | 命令行用法 — 完整参数表、输出结构、报告解读 |
| [docs/api.md](docs/api.md) | Python API — 模块接口、数据类型、工具函数 |
| [docs/algorithm.md](docs/algorithm.md) | 算法说明 — 七步管线详解、调试可视化 |
| [docs/formats.md](docs/formats.md) | 文件格式 — ARC/XYZ/MOL/SDF 详解、扩展新格式 |
| [docs/examples.md](docs/examples.md) | 示例 — 典型分子转换案例、批量处理脚本 |

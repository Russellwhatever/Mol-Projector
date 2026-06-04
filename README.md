# mol-projector

[![中文](https://img.shields.io/badge/README-中文-red)](README.zh-CN.md)

Convert 3D molecular structures to 2D MOL format with automatic bond order determination.

Given a molecule with only atom positions (3D coordinates), `mol-projector` determines single, double, and triple bonds using valence rules, handles special cases (carbenes, conjugated systems, formal charges), and outputs MOL files with correct bond information.

## Install

```bash
# From source
cd mol-projector
pip install -e .

# With conda
conda env create -f env.yaml
conda activate mol-projector
```

详见 [docs/install.md](docs/install.md)。

## Usage

```bash
# Convert all .arc files in a directory
mol-projector -i /path/to/arc/files -o /path/to/output

# With debug visualization
mol-projector -i . -o ./converted -d ./debug_vis

# Specify input format
mol-projector -i . -o ./converted --format xyz
```

详见 [docs/usage.md](docs/usage.md)（完整参数表、输出结构、报告解读）。

## Python API

```python
from mol_projector.io import read_structure
from mol_projector.bond import get_bond_order

# Read a structure file
struc_list = read_structure('molecule.arc', fmt='arc')
struc = struc_list[0]

# Determine bond orders
status, rdmol, smiles = get_bond_order(struc, output_smiles=True)
# status: True (normal), '_e' (electron-modified), '_cb' (carbene),
#         '_cs' (conjugated system), False (failed)

# Write output
from mol_projector.io.mol_io import write_mol_with_fig
write_mol_with_fig(rdmol, 'output.mol', 'output.png')
```

详见 [docs/api.md](docs/api.md)（完整模块接口、数据类型说明）。

## Algorithm

详见 [docs/algorithm.md](docs/algorithm.md)（七步管线详解、调试可视化）。

The bond-order determination pipeline:

1. **Initial bond matrix** — computed from interatomic distances using element-pair cutoff radii (`COMBINE2CUTOFF`)
2. **Multi-bond generation** — upgrade single bonds to double/triple based on valence requirements
3. **Spanning tree** — break cycles to establish base connectivity
4. **Ring closure** — re-add bonds for 5-7 membered and 3-4 membered rings
5. **Multi-bond (second pass)** — final double/triple bond assignment
6. **Electron modifications** — handle formal charges, carbenes (diradicals), and conjugated pi-systems

## Supported Formats

| Format | Extension | Status |
|--------|-----------|--------|
| ARC (BIOVIA/LASP) | `.arc` | Supported |
| XYZ | `.xyz` | Supported |
| MOL (MDL Molfile) | `.mol` | Supported |
| SDF | `.sdf` | Supported |

详见 [docs/formats.md](docs/formats.md)（格式详解、示例、如何扩展新格式）。

## Tests

```bash
pip install pytest
pytest tests/ -v
```

## Verified Examples

| Molecule | Input | SMILES | Bond Orders |
|----------|-------|--------|-------------|
| Water (H2O) | ARC | `[H]O[H]` | O-H single bonds |
| Formaldehyde (H2CO) | ARC | `[H]C([H])=O` | C=O double bond |

详见 [docs/examples.md](docs/examples.md)（批量转换、片段拆分示例）。

## Fragment Decomposition

```bash
# Decompose molecule into ring-based fragments
mol-fragment molecule.arc -f arc -o ./report
```

详见 [docs/examples.md](docs/examples.md) 片段拆分章节。

---

## 文档目录

| 文档 | 说明 |
|------|------|
| [docs/install.md](docs/install.md) | 安装指南 — pip / conda 安装、环境配置、验证 |
| [docs/usage.md](docs/usage.md) | 命令行用法 — 完整参数表、输出结构、报告解读 |
| [docs/api.md](docs/api.md) | Python API — 模块接口、数据类型、工具函数 |
| [docs/algorithm.md](docs/algorithm.md) | 算法说明 — 七步管线详解、调试可视化 |
| [docs/formats.md](docs/formats.md) | 文件格式 — ARC/XYZ/MOL/SDF 详解、扩展新格式 |
| [docs/examples.md](docs/examples.md) | 示例 — 典型分子转换、片段拆分、批量处理脚本 |

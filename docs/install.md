# 安装指南

## 环境要求

- Python >= 3.10

## 依赖

| 包 | 用途 | 是否必需 |
|---|---|---|
| `numpy>=1.21` | 矩阵运算、坐标处理 | 是 |
| `rdkit>=2023.03` | MOL 输出、SMILES 生成 | 是 |
| `tqdm>=4.60` | 进度条 | 是 |
| `networkx>=2.8` | 结构可视化 | 否（`[viz]` 可选） |
| `matplotlib>=3.5` | 结构可视化 | 否（`[viz]` 可选） |

## pip 安装

```bash
cd mol-projector

# 最小安装（不含可视化模块）
pip install -e .

# 完整安装（含可视化模块）
pip install -e ".[viz]"

# 安装开发依赖
pip install -e ".[viz]" pytest
```

## conda 安装

```bash
# 从 env.yaml 创建环境
conda env create -f env.yaml

# 激活环境
conda activate mol-projector

# 在环境中以开发模式安装
pip install -e .
```

## 验证安装

```bash
python -c "from mol_projector.io import read_structure; print('安装成功')"
```

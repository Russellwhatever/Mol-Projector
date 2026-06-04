# mol-projector 文档

## 目录

| 文档 | 内容 |
|------|------|
| [安装指南](install.md) | pip / conda 安装与环境配置 |
| [命令行用法](usage.md) | CLI 参数说明、使用示例、输出结构 |
| [Python API](api.md) | 模块接口与代码示例 |
| [算法说明](algorithm.md) | 键级判定管线的详细步骤 |
| [文件格式](formats.md) | 支持的输入格式及如何扩展 |
| [示例](examples.md) | 典型分子转换案例 |

## 项目概览

`mol-projector` 将仅含三维坐标的分子结构转换为带有正确键级（单键/双键/三键）信息的 MOL 文件。核心功能是从裸坐标推断化学键的键级。

项目结构：

```
mol-projector/
├── src/mol_projector/
│   ├── core/       # 数据模型：Atom（原子）, Struc（结构）, 元素周期表
│   ├── bond/       # 键级判定核心算法（常量、化合价、距离、连通性、电子修正、管线）
│   ├── io/         # 文件读写：ARC, XYZ, MOL, SDF
│   ├── viz/        # 可选：结构可视化（networkx + matplotlib）
│   └── cli/        # 命令行入口
├── tests/          # 单元测试（25 项）
├── docs/           # 文档
└── README.md       # 项目说明（英文 / 中文）
```

# 命令行用法

## 基本语法

```
mol-projector [-i 输入目录] [-o 输出目录] [-d 调试目录] [-f 格式] [-h]
```

## 参数说明

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--input` | `-i` | `.`（当前目录） | 输入目录路径，递归扫描 |
| `--output` | `-o` | `./converted` | 输出根目录 |
| `--debug` | `-d` | 不启用 | 调试可视化输出目录。只写 `-d` 则默认 `./debug_vis` |
| `--format` | `-f` | `arc` | 输入文件格式：`arc`、`xyz`、`mol`、`sdf` |
| `--help` | `-h` | — | 显示帮助信息 |

## 使用示例

```bash
# 转换当前目录下所有 .arc 文件
mol-projector

# 指定输入和输出目录
mol-projector -i /data/arc_files -o ./results

# 开启调试可视化，输出每一步的键连矩阵变化图
mol-projector -i . -o ./converted -d

# 自定义调试输出目录
mol-projector -i . -o ./converted -d ./my_debug

# 转换 XYZ 格式文件
mol-projector -f xyz -i /path/to/xyz/files -o ./xyz_output
```

## 输出目录结构

```
converted/
├── success/          # 成功
│   ├── mol/          #   转换后的 .mol 文件（含键级信息）
│   ├── fig/          #   2D 结构图 .png
│   └── smiles/       #   SMILES 字符串 .smiles
├── fail/             # 键级判定失败：保留原文件副本
├── error/            # 解析出错：保留原文件副本
└── report.txt        # 转换报告
```

## 转换报告示例

```
Conversion Report
================
Total files: 100
Successful conversions: 92
Failed files: 5
Error files: 3
Success rate: 92.00%

Detailed statistics:
- Electron modified (_e): 3       ← 电子修正
- Carbene modified (_cb): 1       ← 卡宾修正
- Conjugated system modified (_cs): 2  ← 共轭体系修正
- Standard conversion: 86

Output directories:
- Success files: ./converted/success
- Failed files: ./converted/fail
- Error files: ./converted/error
```

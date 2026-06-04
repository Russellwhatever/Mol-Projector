# 支持的文件格式

## ARC（BIOVIA / LASP 存档格式）

**扩展名**：`.arc`

**格式说明**：LASP 使用的 BIOVIA 存档格式，也是本项目的原生格式。每个结构以 `Energy` 行开始，包含原子坐标块（`CORE` 或 `XXXX` 行）、可选的 `PBC` 晶格参数行和 `!DATE` 标签行。

**示例**：

```
 Energy -40.1234 0.0012
 O      0.000000000    0.000000000    0.117000000  CORE  0  0  0  0  0
 H      0.000000000    0.757000000   -0.469000000  CORE  0  0  0  0  0
 H      0.000000000   -0.757000000   -0.469000000  CORE  0  0  0  0  0
 end
```

**解析方式**：本项目自带纯 Python 解析器（`io/arc_parser.py`），不依赖 laspy。同时保留了与 laspy 的接口：可通过上层调用 `laspy.basic.allstr.AllStr.read_arc()` 读取后手动构造 `Struc` 对象传入管线。

## XYZ

**扩展名**：`.xyz`

**格式说明**：标准 XYZ 格式。第一行为原子数，第二行为注释行，后续每行为元素符号和 x y z 坐标。

**示例**：

```
3
water molecule, energy=-40.12
O    0.000000    0.000000    0.117000
H    0.000000    0.757000   -0.469000
H    0.000000   -0.757000   -0.469000
```

**扩展 XYZ**：支持 ASE 风格扩展格式，可从注释行解析 `energy=` 和 `Lattice=` 字段。

## MOL（MDL Molfile）

**扩展名**：`.mol`

**格式说明**：MDL Molfile V2000 格式。通过 RDKit `Chem.MolFromMolFile` 读取，再转换为内部 `Struc` 对象。

**输入**：带有键级信息的 MOL 文件，键级作为参考值保留。

**输出**：管线输出的 MOL 文件包含正确的键级、形式电荷和自由基信息。

## SDF

**扩展名**：`.sdf`

**格式说明**：MDL SDF 格式，单文件可包含多个分子。通过 RDKit `Chem.SDMolSupplier` 逐分子读取，每个转换为独立的 `Struc` 对象。

## 扩展新格式

通过注册器装饰器即可添加新格式解析器：

```python
from mol_projector.io import register_reader
from mol_projector.core import Struc

@register_reader('pdb')
def read_pdb(path):
    """PDB 格式解析器"""
    struc_list = []
    # ... 解析逻辑 ...
    return struc_list
```

注册后 CLI 和 API 均可直接使用：

```bash
mol-projector -f pdb -i /path/to/pdb/files
```

```python
from mol_projector.io import read_structure
struc = read_structure('protein.pdb', fmt='pdb')
```

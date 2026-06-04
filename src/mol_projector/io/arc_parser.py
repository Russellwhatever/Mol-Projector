"""ARC format parser (BIOVIA/LASP archive format)."""

import os
import re
import numpy as np
from . import register_reader
from ..core import Struc
from ..core.periodictable import ele_dict


@register_reader('arc')
def read_arc(filename: str):
    """Parse a BIOVIA/LASP ARC file.

    Args:
        filename: Path to the .arc file

    Returns:
        List of Struc objects
    """
    result = []
    if not os.path.isfile(filename):
        print(f'--No {filename} at {os.getcwd()}--')
        return result

    with open(filename) as f:
        currentStr = -1
        lastline = ""
        second_lastline = ""

        for line in f:
            if 'Energy' in line or 'React' in line:
                result.append(Struc())
                currentStr += 1

                parts = line.split()
                if parts == ["Energy"]:
                    result[currentStr].energy = 0.0
                elif len(parts) >= 3:
                    try:
                        result[currentStr].energy = float(parts[1])
                    except ValueError:
                        result[currentStr].energy = float(parts[2]) if len(parts) > 2 else 0.0
                elif '****' in line:
                    result[currentStr].energy = 100.0
                else:
                    result[currentStr].energy = 0.0
                try:
                    result[currentStr].max_f = float(parts[-1]) if len(parts) >= 3 else 0.0
                except (ValueError, IndexError):
                    result[currentStr].max_f = 0

            elif 'CORE' in line or 'XXXX' in line:
                li = line.split()
                ele = re.sub(r"\d+", "", li[0])
                ele = ele_dict[ele]
                if len(li) == 10:
                    xyz = [float(x) for x in li[1:4]]
                elif len(li) < 10:
                    xyz = [float(line[5:20]), float(line[20:35]), float(line[35:50])]
                result[currentStr].add_atom(ele, xyz)

            elif 'PBC' in line and 'ON' not in line:
                result[currentStr].abc = np.array(
                    [float(x) for x in line.split()[1:7]])

            elif '!DATE' in line:
                result[currentStr].label = line[5:-1]
                if len(line.split()) >= 2:
                    result[currentStr].molecule_name = line.split()[-1]
                    result[currentStr].label = ' '.join(line.split()[1:])

            second_lastline = lastline
            lastline = line

        if lastline.strip():
            if 'end' not in lastline and "PBC=ON" not in lastline:
                if len(result) > 0:
                    result.pop(-1)
        else:
            lastline = second_lastline
            if 'end' not in lastline and "PBC=ON" not in lastline:
                if len(result) > 0:
                    result.pop(-1)

    return result

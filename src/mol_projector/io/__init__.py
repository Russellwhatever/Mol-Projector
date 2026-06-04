"""I/O module: format readers and writers."""

from typing import Dict, Callable, List

_READERS: Dict[str, Callable] = {}

def register_reader(fmt: str):
    """Decorator to register a format reader."""
    def decorator(func):
        _READERS[fmt.lower()] = func
        return func
    return decorator

def read_structure(path: str, fmt: str = 'arc') -> List:
    """Read molecular structure(s) from a file.

    Args:
        path: Path to the input file
        fmt: Format name ('arc', 'xyz', 'mol', 'sdf')

    Returns:
        List of Struc objects
    """
    reader = _READERS.get(fmt.lower())
    if reader is None:
        raise ValueError(f"Unknown format: {fmt}. Supported: {list(_READERS.keys())}")
    return reader(path)

def list_supported_formats() -> List[str]:
    """Return list of supported input formats."""
    return list(_READERS.keys())

# Import readers to trigger registration
from . import arc_parser
from . import xyz_parser
from . import mol_io

from .readpath import ReadPath

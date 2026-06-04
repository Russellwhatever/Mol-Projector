"""Directory scanner for molecule files."""

import os


class ReadPath(dict):
    """Scan a directory (recursively) for files of a given format.

    Maps molecule names to file paths.

    Usage:
        name2path = ReadPath('/path/to/dir', format='arc')
        for name, path in name2path.items():
            ...
    """

    def __init__(self, input_path: str, fmt: str = 'arc', get_name_f=None):
        super().__init__()
        self.format = fmt
        self.get_name_f = get_name_f
        self._read_path(input_path)

    def _read_path(self, input_path):
        if isinstance(input_path, str):
            if os.path.isfile(input_path) and input_path.endswith(f'.{self.format}'):
                self.update({self.file2name(input_path): input_path})
            elif os.path.isdir(input_path):
                for file_name in os.listdir(input_path):
                    self._read_path(os.path.join(input_path, file_name))
        elif isinstance(input_path, list):
            for p in input_path:
                self._read_path(p)

    def file2name(self, file_path: str) -> str:
        if self.get_name_f:
            return self.get_name_f(file_path)
        else:
            return file_path.split('/')[-1][:-(len(self.format) + 1)]

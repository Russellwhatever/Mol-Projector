#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Molecular structure conversion script.

Converts 3D molecular structure files to MOL format with bond order determination.

Usage:
    mol-projector [options]

Options:
    -i, --input DIR         Input directory path (default: current directory)
    -o, --output DIR        Output directory path (default: ./converted)
    -d, --debug DIR         Debug visualization directory path (default: None)
    -f, --format FMT        Input file format (default: arc)
    -h, --help              Show help message

Examples:
    mol-projector
    mol-projector -i /path/to/arc/files -o /path/to/output
    mol-projector --input . --output ./results --debug ./debug
    mol-projector -f xyz -i /path/to/xyz/files
"""

import argparse
import os
import sys
import shutil
from tqdm import tqdm

from mol_projector.io import ReadPath, read_structure, list_supported_formats
from mol_projector.io.mol_io import write_mol_with_fig
from mol_projector.bond.pipeline import get_bond_order
from mol_projector.core import Struc


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert 3D molecule to 2D MOL format with bond order determination",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '-i', '--input',
        type=str,
        default='.',
        help='Input directory (default: .)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='./converted',
        help='Output directory (default: ./converted)'
    )

    parser.add_argument(
        '-d', '--debug',
        type=str,
        nargs='?',
        const='./debug_vis',
        default=None,
        help='Debug visualization directory (default: ./debug_vis when specified, None when not)'
    )

    parser.add_argument(
        '-f', '--format',
        type=str,
        default='arc',
        choices=list_supported_formats(),
        help=f'Input file format (default: arc). Supported: {list_supported_formats()}'
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    input_dir = os.path.abspath(args.input)
    converted_dir = os.path.abspath(args.output)
    debug_dir = args.debug
    fmt = args.format

    if not os.path.exists(input_dir):
        print(f"Input directory does not exist: {input_dir}")
        sys.exit(1)

    success_dir = os.path.join(converted_dir, 'success')
    fail_dir = os.path.join(converted_dir, 'fail')
    error_dir = os.path.join(converted_dir, 'error')

    os.makedirs(success_dir, exist_ok=True)
    os.makedirs(fail_dir, exist_ok=True)
    os.makedirs(error_dir, exist_ok=True)
    if debug_dir:
        os.makedirs(debug_dir, exist_ok=True)

    # Clear output directories
    dirs_to_clear = [success_dir, fail_dir, error_dir]
    if debug_dir:
        dirs_to_clear.append(debug_dir)

    for dir_path in dirs_to_clear:
        for file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, file)
            if os.path.isfile(file_path):
                os.remove(file_path)

    # Read input files
    try:
        name2path = ReadPath(input_dir, fmt=fmt)
    except Exception as e:
        print(f"Error: failed to read {fmt} files: {e}")
        sys.exit(1)

    if not name2path:
        print(f"WARNING: no .{fmt} files found at {input_dir}")
        sys.exit(0)

    # Statistics
    success_num = 0
    total_num = len(name2path)
    e_num = 0
    cb_num = 0
    cs_num = 0
    error_num = 0

    for name, file_path in tqdm(name2path.items(), desc=f"Converting {fmt} files"):
        try:
            struc_list = read_structure(file_path, fmt=fmt)
            if not struc_list:
                save_path = os.path.join(error_dir, f"{name}.{fmt}")
                shutil.copy(file_path, save_path)
                error_num += 1
                continue

            struc = struc_list[0]

            # Debug visualization
            debug_file_path = None
            if debug_dir:
                debug_file_path = os.path.join(debug_dir, f"{name}.png")

            # Run bond-order determination
            status, mol, smiles = get_bond_order(
                struc,
                output_smiles=True,
                visualize_path=debug_file_path
            )

            if status:
                # Write output
                mol_dir = os.path.join(success_dir, 'mol')
                fig_dir = os.path.join(success_dir, 'fig')
                smiles_dir = os.path.join(success_dir, 'smiles')

                if isinstance(status, str):
                    name_with_suffix = name + status
                else:
                    name_with_suffix = name

                mol_path = os.path.join(mol_dir, f"{name_with_suffix}.mol")
                fig_path = os.path.join(fig_dir, f"{name_with_suffix}.png")
                write_mol_with_fig(mol, mol_path, fig_path)

                if smiles:
                    os.makedirs(smiles_dir, exist_ok=True)
                    smiles_path = os.path.join(smiles_dir, f"{name_with_suffix}.smiles")
                    with open(smiles_path, 'w') as f:
                        f.write(smiles)

                success_num += 1
                if status == '_e':
                    e_num += 1
                elif status == '_cb':
                    cb_num += 1
                elif status == '_cs':
                    cs_num += 1
            else:
                save_path = os.path.join(fail_dir, f"{name}.{fmt}")
                shutil.copy(file_path, save_path)

        except Exception as e:
            print(f"Error at {file_path}: {str(e)}")
            save_path = os.path.join(error_dir, f"{name}.{fmt}")
            try:
                shutil.copy(file_path, save_path)
            except Exception:
                pass
            error_num += 1

    # Generate report
    if total_num > 0:
        success_rate = success_num / total_num
    else:
        success_rate = 0

    report_content = f"""Conversion Report
================

Total files: {total_num}
Successful conversions: {success_num}
Failed files: {total_num - success_num - error_num}
Error files: {error_num}
Success rate: {success_rate:.2%}

Detailed statistics:
- Electron modified (_e): {e_num}
- Carbene modified (_cb): {cb_num}
- Conjugated system modified (_cs): {cs_num}
- Standard conversion: {success_num - e_num - cb_num - cs_num}

Output directories:
- Success files: {success_dir}
- Failed files: {fail_dir}
- Error files: {error_dir}"""

    if debug_dir:
        report_content += f"\n- Visualized structures at: {debug_dir}"

    report_content += "\n"

    report_path = os.path.join(converted_dir, 'report.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report_content)

    print(f"\n{success_rate:.2%} of {total_num} successfully converted")
    print(f"e_modified {e_num}, carbene {cb_num}, converge_system {cs_num}, normal {success_num - e_num - cb_num - cs_num}")
    print(f"Detailed report at: {report_path}")


if __name__ == '__main__':
    main()

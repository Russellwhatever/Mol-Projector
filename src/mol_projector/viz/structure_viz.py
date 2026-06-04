"""Molecular structure visualization using networkx and matplotlib."""

import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


def _get_next_filename(base_path):
    """Auto-increment filename to avoid overwriting."""
    if not os.path.exists(base_path):
        return base_path

    directory = os.path.dirname(base_path)
    filename = os.path.basename(base_path)
    name, ext = os.path.splitext(filename)

    if '_' in name:
        base_name, num = name.rsplit('_', 1)
        if num.isdigit():
            name = base_name
            current_num = int(num)
        else:
            current_num = 0
    else:
        current_num = 0

    while True:
        current_num += 1
        new_filename = f"{name}_{current_num}{ext}"
        new_path = os.path.join(directory, new_filename)
        if not os.path.exists(new_path):
            return new_path


def visualize_structure(struc, save_path='/tmp/mol_projector_structure.png'):
    """Visualize molecular structure with atom indices, element labels, and bond orders.

    Args:
        struc: Struc object with bond_matrix
        save_path: Path to save the PNG image
    """
    save_path = _get_next_filename(save_path)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    adj_matrix = (struc.bond_matrix != 0).astype(int)
    bond_matrix = struc.bond_matrix

    G = nx.Graph()

    for i, atom in enumerate(struc.atom):
        G.add_node(i, element=atom.ele_symb)

    edge_labels = {}
    for i in range(len(bond_matrix)):
        for j in range(i + 1, len(bond_matrix)):
            if adj_matrix[i, j] > 0:
                G.add_edge(i, j)
                bond_order = int(bond_matrix[i, j])
                edge_labels[(i, j)] = str(bond_order)

    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color='lightblue',
                           node_size=1000, alpha=0.6)

    edge_colors = []
    edge_widths = []
    for (i, j) in G.edges():
        bond_order = int(bond_matrix[i, j])
        if bond_order == 1:
            edge_colors.append('gray')
            edge_widths.append(1)
        elif bond_order == 2:
            edge_colors.append('blue')
            edge_widths.append(2)
        elif bond_order == 3:
            edge_colors.append('red')
            edge_widths.append(3)
        else:
            edge_colors.append('black')
            edge_widths.append(4)

    nx.draw_networkx_edges(G, pos, edge_color=edge_colors,
                           width=edge_widths, alpha=0.6)

    labels = {i: f"{i}\n{atom.ele_symb}" for i, atom in enumerate(struc.atom)}
    nx.draw_networkx_labels(G, pos, labels, font_size=10)

    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,
                                 font_size=8, font_color='black')

    plt.title(f"Molecular Structure Visualization\nTotal atoms: {struc.natom}")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Structure visualization saved to {save_path}")

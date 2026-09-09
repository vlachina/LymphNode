import tyssue
import os
import pickle
from tyssue import History
from tyssue.draw import sheet_view
from tyssue.draw.plt_draw import quick_edge_draw
from tyssue.io import obj
import matplotlib.pyplot as plt
import pathlib
import random
import numpy as np
from scipy.stats import skewnorm#
from matplotlib import cm
from matplotlib.cm import ScalarMappable
import logging
import pandas as pd
from scipy.stats import truncnorm
from tyssue.draw.plt_draw import draw_edge, draw_vert
from matplotlib.colors import ListedColormap


def identify_boundary_layers(cellmap, max_layers):

    """
    Identifies which faces/edges/vertices lie within a given number of layers from the tissue boundary,
    so boundary effects can be avoided or a custom boundary condition applied.
    """

    # Keeping track of the boundary faces, inside faces, boundary edges, inside edges, inside vertices
    boundary_faces = []
    inside_faces = []
    boundary_edges = []
    inside_edges = []
    inside_vertices = []

    # Keeping track of outside vertices
    outside_vertices = set()

    # Keeping a working DataFrame for included faces and edges
    working_faces_df = cellmap.face_df.copy()
    working_edges_df = cellmap.edge_df.copy()

    for layer in range(max_layers):
        current_edges_without_opposites = []

        # Identifying edges without opposites in the current working set
        for edge_idx, edge_row in working_edges_df.iterrows():
            if edge_row['face'] in working_faces_df.index:
                srce, trgt = edge_row['srce'], edge_row['trgt']
                opposite_edges = working_edges_df[
                    (working_edges_df['srce'] == trgt) &
                    (working_edges_df['trgt'] == srce)
                    ]
                if opposite_edges.empty:
                    current_edges_without_opposites.append(edge_idx)

        # Excluding faces associated with these edges
        faces_to_exclude = working_edges_df.loc[current_edges_without_opposites, 'face'].unique()

        # Excluding other edges associated with the given faces
        other_edges_to_exclude = working_edges_df[working_edges_df['face'].isin(faces_to_exclude)].index.tolist()

        # Adding these faces and edges to the boundary list
        boundary_faces.extend(faces_to_exclude)
        boundary_edges.extend(other_edges_to_exclude)  # Exclude all edges related to these faces

        # Excluding these faces from the working DataFrame
        working_faces_df = working_faces_df[~working_faces_df.index.isin(faces_to_exclude)]
        working_edges_df = working_edges_df[~working_edges_df['face'].isin(faces_to_exclude)]

        # If no more faces can be excluded, stop the loop
        if working_faces_df.empty:
            break

    # Identifying inside edges and faces
    inside_faces = working_faces_df.index.tolist()  # The remaining faces
    inside_edges = working_edges_df.index.tolist()  # The remaining edges

    # Identifying inside vertices by looping over the inside edges
    for edge_id in inside_edges:
        edge = cellmap.edge_df.loc[edge_id]
        srce_vertex = edge['srce']
        trgt_vertex = edge['trgt']

        # Adding the vertices to the inside vertices list
        inside_vertices.append(srce_vertex)
        inside_vertices.append(trgt_vertex)

    # Removing duplicates in the inside vertices list
    inside_vertices = list(set(inside_vertices))

    # Identifying outside vertices (those not in inside vertices)
    for edge_id in boundary_edges:
        edge = cellmap.edge_df.loc[edge_id]
        srce_vertex = edge['srce']
        trgt_vertex = edge['trgt']

        # Adding boundary vertices to the outside vertices list if they are not inside vertices
        if srce_vertex not in inside_vertices:
            outside_vertices.add(srce_vertex)
        if trgt_vertex not in inside_vertices:
            outside_vertices.add(trgt_vertex)

    # Converting the set of outside vertices back to a list
    outside_vertices = list(outside_vertices)

    # Identifying outside edges (edges touching outside vertices)
    outside_edges = []
    for edge_id, edge in cellmap.edge_df.iterrows():
        if edge['srce'] in outside_vertices or edge['trgt'] in outside_vertices:
            outside_edges.append(edge_id)

    return boundary_edges, boundary_faces, inside_edges, outside_edges, inside_faces, inside_vertices, outside_vertices


def view(cellmap, geom, show_axes=True, xlabel="x", ylabel="y", xlim=None, ylim=None):
    """
    Visualises the cell map with optional axes and labels.

    Parameters:
        show_axes (bool): If True, shows axis ticks and labels.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        xlim (tuple): Tuple specifying (xmin, xmax) for x-axis limits.
        ylim (tuple): Tuple specifying (ymin, ymax) for y-axis limits.
    """
    geom.update_all(cellmap)

    ecolor = "blue"

    fig, ax = sheet_view(
        cellmap,
        mode="2D",
        face={"visible": False},
        edge={"head_width": 0.0, "color": ecolor, "width": 2, "alpha": 1.0},
        vert={"visible": True, "s": 10, "color": "red", "alpha": 1.0}
    )
    fig.set_size_inches(15, 15)

    if show_axes:
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    # Set axis limits if provided
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    return fig, ax


def highlight_vertices(cellmap, geom, chosen_vert_ids, default_color='red', highlight_color='lime', default_size=10, highlight_size=40, show_axes=True):
    
    """
    Plots the tissue with all edges shown, highlighting specific vertices
    (e.g. a vertex just created by division or a collapsed edge's merge point)
    in a different colour and size from the rest.

    - chosen_vert_ids: Vertex ID(s) to highlight.
    - default_color, default_size: Appearance of all non-highlighted vertices.
    - highlight_color, highlight_size: Appearance of the highlighted vertices.
    - show_axes: If True, shows x/y axis labels; if False, hides tick marks entirely.
    """

    # Updating geometry before plotting
    geom.update_all(cellmap)

    # Setting default vertex appearance
    cellmap.vert_df['color'] = default_color
    cellmap.vert_df['size'] = default_size

    # Highlighting specified vertices
    for v_id in chosen_vert_ids:
        if v_id in cellmap.vert_df.index:
            cellmap.vert_df.loc[v_id, 'color'] = highlight_color
            cellmap.vert_df.loc[v_id, 'size'] = highlight_size


    # Plotting with consistent visual style
    fig, ax = sheet_view(
        cellmap,
        mode="2D",
        face={"visible": False},
        edge={"head_width": 0.0, "color": "blue", "width": 2, "alpha": 1.0},
        vert={
            "visible": True,
            "s": cellmap.vert_df['size'],
            "color": cellmap.vert_df['color'],
            "alpha": 1.0
        }
    )

    fig.set_size_inches(15, 15)
    if show_axes:
        ax.set_xlabel("x", fontsize=14)
        ax.set_ylabel("y", fontsize=14)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    return fig, ax



def load_simulation_state(folder_path, filename='final_cellmap_state.pkl'):

    """
    Loads a saved pkl file 
    """

    file_path = os.path.join(folder_path, filename)

    # Loading the cellmap object using pickle
    with open(file_path, 'rb') as f:
        cellmap = pickle.load(f)

    print(f"Simulation state loaded from {file_path}")
    return cellmap


def save_simulation_state(cellmap, folder_path, filename='final_cellmap_state.pkl'):

    """
    Saves current model state as a pkl file 
    """

    # Ensuring the folder exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Constructing the full file path
    file_path = os.path.join(folder_path, filename)

    # Saving the cellmap object using pickle
    with open(file_path, 'wb') as f:
        pickle.dump(cellmap, f)

    print(f"Simulation state saved at {file_path}")


def highlight_edge_on_cellmap(cellmap_init, edge_id=None, figsize=(15, 15),
                              vert_col='orange', edge_highlight_col='red',
                              base_edge_col='black', base_vert_col='black', vert_size=50,
                              save_path=None, show_axes=False,
                              xlim=None, ylim=None,
                              plot_vertices=True, show_figure=True):
    """
    Plots a Tyssue cellmap with optional edge highlighting.

    Parameters:
    - cellmap_init: The original cellmap to copy and draw.
    - edge_id: (Optional) The edge index to highlight.
    - figsize: Tuple defining the figure size.
    - vert_col: Colour of the top-layer vertices.
    - edge_highlight_col: Colour for the highlighted edge (used in custom colormap).
    - base_edge_col: Colour of the underlying edges (for outlines).
    - base_vert_col: Colour of the underlying vertices (for outlines).
    - save_path: Optional path to save the figure.
    - plot_vertices: If True, plot vertices; if False, skip vertex plotting.
    - show_figure: If True, call plt.show(); if False, don't display the figure.
    """

   

    # Copy the cellmap
    cellmap = cellmap_init.copy()
    coords = cellmap.coords[:2]  # 2D coordinates

    # Set edge colour map indices and widths
    cellmap.edge_df["color_edges"] = 0
    cellmap.edge_df["edge_width"] = 3

    if edge_id is not None:
        if edge_id in cellmap.edge_df.index:
            cellmap.edge_df.loc[edge_id, "color_edges"] = 1
            cellmap.edge_df.loc[edge_id, "edge_width"] = 4
        else:
            print(f"Warning: Edge ID {edge_id} not found in edge_df")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Draw base (thicker black) edges
    draw_edge(
        cellmap,
        coords=coords,
        ax=ax,
        color=base_edge_col,
        width=cellmap.edge_df["edge_width"] + 2,
        head_width=0.0,
        alpha=1.0,
        zorder=2
    )

    # Draw all edges uniformly in yellow
    draw_edge(
        cellmap,
        coords=coords,
        ax=ax,
        color="yellow",
        width=cellmap.edge_df["edge_width"],
        head_width=0.0,
        alpha=1.0,
        zorder=2
    )

    # Draw the highlighted edge last, directly, so its parallel (opposite-direction)
    # twin can never be drawn on top of it
    if edge_id is not None and edge_id in cellmap.edge_df.index:
        srce, trgt = cellmap.edge_df.loc[edge_id, ["srce", "trgt"]]
        x_vals = [cellmap.vert_df.loc[srce, coords[0]], cellmap.vert_df.loc[trgt, coords[0]]]
        y_vals = [cellmap.vert_df.loc[srce, coords[1]], cellmap.vert_df.loc[trgt, coords[1]]]
        ax.plot(x_vals, y_vals, color=edge_highlight_col,
                linewidth=cellmap.edge_df.loc[edge_id, "edge_width"], zorder=3)

    # Draw vertices (only if plot_vertices is True)
    if plot_vertices:
        # Draw vertex outlines (underlay)
        draw_vert(
            cellmap,
            coords=coords,
            ax=ax,
            color=base_vert_col,
            s=vert_size,
            alpha=1.0,
            zorder=3
        )

        # Draw foreground vertices
        draw_vert(
            cellmap,
            coords=coords,
            ax=ax,
            color=vert_col,
            s=vert_size - 40,
            alpha=1.0,  
            zorder=4
        )

    # Format
    if show_axes:
        ax.set_xlabel("X", fontsize=14)
        ax.set_ylabel("Y", fontsize=14)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    ax.set_aspect('equal')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    
    if show_figure:
        plt.show()
    elif not show_figure and save_path is None:
        plt.close(fig)  # Close the figure if not showing and not saving
    
    return fig, ax
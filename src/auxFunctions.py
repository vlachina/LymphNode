import tyssue
import os
import pickle
from tyssue import History
from tyssue.draw import sheet_view
from tyssue.draw.plt_draw import quick_edge_draw
from tyssue.io import obj
import matplotlib.pylab as plt
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


def sample_truncated_normal(mean, std_dev, lower_bound, upper_bound, size=1):
    a, b = (lower_bound - mean) / std_dev, (upper_bound - mean) / std_dev
    return truncnorm.rvs(a, b, loc=mean, scale=std_dev, size=size)


def length_elasticity_range(cellmap, lower_FRC_le, higher_FRC_le, lower_ECM_le, higher_ECM_le):
    # Create new columns if not already present
    cellmap.edge_df['length_elasticity_ECM'] = cellmap.edge_df['length_elasticity']
    cellmap.edge_df['length_elasticity_FRC'] = cellmap.edge_df['length_elasticity']

    for edge in range(len(cellmap.edge_df)):
        newFRCValue = random.uniform(lower_FRC_le, higher_FRC_le)
        newECMValue = random.uniform(lower_ECM_le, higher_ECM_le)

        # Update edge-specific parameters
        cellmap.edge_df.at[edge, 'length_elasticity_ECM'] = newECMValue
        cellmap.edge_df.at[edge, 'length_elasticity_FRC'] = newFRCValue
        cellmap.edge_df.at[edge, 'length_elasticity'] = newFRCValue + newECMValue

    return cellmap


def length_elasticity_ECM_range(cellmap, lower_ECM_le, higher_ECM_le):
    """
    Assigns a random ECM elasticity value (between lower_ECM_le and higher_ECM_le)
    to each edge and stores it in 'length_elasticity_ECM'.
    """

    for edge in cellmap.edge_df.index:
        newECMValue = random.uniform(lower_ECM_le, higher_ECM_le)
        cellmap.edge_df.loc[edge, 'length_elasticity'] = newECMValue

    return cellmap

def line_tension_FRC_range(cellmap, lower_FRC_le, higher_FRC_le):
    """
    Assigns a random FRC elasticity value (between lower_FRC_le and higher_FRC_le)
    to each edge and stores it in 'length_elasticity_FRC'.
    """

    for edge in cellmap.edge_df.index:
        newFRCValue = random.uniform(lower_FRC_le, higher_FRC_le)
        cellmap.edge_df.loc[edge, 'line_tension'] = newFRCValue

    return cellmap


def set_opposite_edges_mechanics(cellmap):

    for edge in range(len(cellmap.edge_df)):
    
    ### find the opposite edge of each edge in edges_list_without_opposites
    
        opposite_edge_indices = cellmap.edge_df[
            (cellmap.edge_df["srce"] == cellmap.edge_df.loc[edge, "trgt"]) &
            (cellmap.edge_df["trgt"] == cellmap.edge_df.loc[edge, "srce"])
        ].index
    
    ### set opposite edges parameters

        if len(opposite_edge_indices):
            # Update the values for the opposite edges
                cellmap.edge_df['length_elasticity_ECM'][opposite_edge_indices] = cellmap.edge_df['length_elasticity_ECM'][edge] 
                cellmap.edge_df['length_elasticity_FRC'][opposite_edge_indices] = cellmap.edge_df['length_elasticity_FRC'][edge] 
                cellmap.edge_df['length_elasticity'][opposite_edge_indices] = cellmap.edge_df['length_elasticity'][edge]

    return cellmap

def line_tension(cellmap, lower_line_tension, higher_line_tension):
    for edge in range(len(cellmap.edge_df)):
        newValue = random.uniform(lower_line_tension, higher_line_tension)
        cellmap.edge_df.loc[edge, 'line_tension'] = newValue 
    return cellmap

def area_elasticity_range(cellmap, lower_area_elasticity, higher_area_elasticity):
    for face in range(len(cellmap.face_df)):
        newValue = random.uniform(lower_area_elasticity, higher_area_elasticity)
        cellmap.face_df.loc[face, 'area_elasticity'] = newValue 
    return cellmap

def prefered_area_range(cellmap, lower_prefered_area, higher_prefered_area):
    for face in range(len(cellmap.face_df)):
        newValue = random.uniform(lower_prefered_area, higher_prefered_area)
        cellmap.face_df.loc[face, 'prefered_area'] = newValue
    return cellmap


def random_island_prefered_area(cellmap, area_large, area_small, probability_large=0.5):
    """
    Assigns preferred area values randomly across the network.

    Parameters:
    - cellmap: The cell map containing FRCs.
    - area_large: The preferred area for large cells.
    - area_small: The preferred area for small cells.
    - probability_large: The probability of assigning a cell the large area.
    """
    for face in range(len(cellmap.face_df)):
        # Randomly assign large or small preferred area based on the specified probability
        if random.random() < probability_large:
            cellmap.face_df.loc[face, 'prefered_area'] = area_large
        else:
            cellmap.face_df.loc[face, 'prefered_area'] = area_small
    return cellmap

def prefered_length_range(cellmap, lower_prefered_length, higher_prefered_length):
    for edge in range(len(cellmap.edge_df)):
        newValue = random.uniform(lower_prefered_length, higher_prefered_length)
        cellmap.edge_df.loc[edge, 'prefered_length'] = newValue 
    return cellmap



def create_frames(
    history,
    output,
    num_frames=None,
    interval=None,
    draw_func=None,
    margin=5,
    **draw_kwds,
):
    """Creates a set of png frames of the recorded history.
   
    Parameters
    ----------
    history : a :class:`tyssue.History` object
    output : path to the output directory
    num_frames : int, the number of frames in the gif
    interval : tuples, define begin and end frame of the gif
    draw_func : a drawing function
         this function must take a `sheet` object as first argument
         and return a `fig, ax` pair. Defaults to quick_edge_draw
         (aka sheet_view with quick mode)
    margin : int, the graph margins in percents, default 5
         if margin is -1, let the draw function decide
    **draw_kwds are passed to the drawing function
    """
    graph_dir = pathlib.Path(output)
    graph_dir.mkdir(parents=True, exist_ok=True)

    x, y = coords = draw_kwds.get("coords", history.sheet.coords[:2])
    sheet0 = history.retrieve(0)
    bounds = sheet0.vert_df[coords].describe().loc[["min", "max"]]
    delta = (bounds.loc["max"] - bounds.loc["min"]).max()
    margin = delta * margin / 100
    xlim = bounds.loc["min", x] - margin, bounds.loc["max", x] + margin
    ylim = bounds.loc["min", y] - margin, bounds.loc["max", y] + margin

    if interval is None:
        start, stop = None, None
    else:
        start, stop = interval[0], interval[1]

    for i, (t, sheet) in enumerate(history.browse(start, stop, num_frames)):
        try:
            fig, ax = sheet_view(sheet, **draw_kwds)
            fig.set_size_inches(20,20)

            if isinstance(ax, plt.Axes) and margin >= 0:
                ax.set(xlim=xlim, ylim=ylim)

            plt.axis('off')
            fig.savefig(graph_dir / f"movie_{i:04d}.png")
        except Exception as e:
            print("Droped frame {i}")
            print(e)

        plt.close()

def exportToMesh(history, dir):
    """Exporting each timepoint to mesh"""
    for i, (t, sheet) in enumerate(history.browse(None, None, None)):
        obj.save_splitted_cells(dir + '/junctions_'+ str(t) +'.obj', sheet, epsilon=0.001)

def draw_vert(cellmap):
    draw_specs = tyssue.config.draw.sheet_spec()
    draw_specs['vert']['visible'] = True
    draw_specs['vert']['color'] = "grey"
    draw_specs['vert']['alpha'] = 0.5
    draw_specs['vert']['s'] = 5
    coords = ['x', 'y']
    fig, ax = sheet_view(cellmap, coords, **draw_specs)
    fig.set_size_inches((10, 10))

    return cellmap


import matplotlib.pyplot as plt


def highlight_vertices(cellmap, geom, chosen_vert_ids, default_color='grey', highlight_color='red', default_size=5, highlight_size=20):
    # Update geometry
    geom.update_all(cellmap)

    # Set draw specs
    draw_specs = tyssue.config.draw.sheet_spec()
    draw_specs['edge']['visible'] = True
    draw_specs['edge']['color'] = 'blue'

    draw_specs['vert']['visible'] = True
    draw_specs['vert']['color'] = default_color
    draw_specs['vert']['alpha'] = 1
    draw_specs['vert']['s'] = default_size

    # Set default colours/sizes
    cellmap.vert_df['color'] = default_color
    cellmap.vert_df['size'] = default_size

    # Highlight
    for chosen_vert_id in chosen_vert_ids:
        if chosen_vert_id in cellmap.vert_df.index:
            cellmap.vert_df.loc[chosen_vert_id, 'color'] = highlight_color
            cellmap.vert_df.loc[chosen_vert_id, 'size'] = highlight_size

    draw_specs['vert']['color'] = cellmap.vert_df['color']
    draw_specs['vert']['s'] = cellmap.vert_df['size']

    fig, ax = sheet_view(cellmap, ['x', 'y'], **draw_specs)
    fig.set_size_inches((10, 10))
    plt.show()



def heatmap_of_edges(cellmap, parameter):

    ## condition should be in the form: cellmap.edge_df["parameter"]

    specs = {
        'face': {
            'visible': False,
        },
        'edge': {
            'visible': True,
            'color': parameter,
            'colormap': 'RdPu',
            'width': 2,
        },
        'vert': {
            'visible': False,
        }
    }

    fig, ax = sheet_view(cellmap, **specs)

    # Create a ScalarMappable with the dynamic range of your data for the colorbar
    sm = ScalarMappable(cmap=cm.RdPu, norm=plt.Normalize(vmin=parameter.min(), vmax=parameter.max()))
    sm.set_array(parameter)  # Set the actual array of values for accurate color mapping

    # Add a colorbar to the plot with the dynamic range
    #cbar = plt.colorbar(sm, ax=ax)

    fig.set_size_inches(20, 20)

    plt.show()  # Display the plot

def colour_face(cellmap, face):
    
    cellmap.face_df["color_faces"] = 0

    cellmap.face_df.loc[face, "color_faces"] = 1

    specs = {
        'face': {
            'visible': True,
            'color': cellmap.face_df['color_faces'],
            'colormap':'PuRd',
        },
        'edge': {
            'visible': False,
        },
        'vert': {
            'visible': False,
        }
    }

    fig, ax = sheet_view(cellmap, **specs)
    fig.set_size_inches((10, 10))

def colour_edge(cellmap, geom, edge):
    """
    Visualises the cell map with solid blue edges and highlights a specific edge in red.

    Parameters:
    - cellmap: Tyssue sheet.
    - geom: Geometry object to update geometry before plotting.
    - edge: The ID(s) of edge(s) to highlight.
    """
    # Update geometry
    geom.update_all(cellmap)

    # Initialise all edge colours and widths
    cellmap.edge_df["color_edges"] = 0  # Default: blue in 'bwr' colormap
    cellmap.edge_df['edge_width'] = 2

    # Highlight selected edge(s)
    cellmap.edge_df.loc[edge, "color_edges"] = 1  # Will become red in 'bwr'
    cellmap.edge_df.loc[edge, 'edge_width'] = 4

    # Prepare specs
    specs = {
        'face': {'visible': False},
        'edge': {
            'visible': True,
            'color': cellmap.edge_df['color_edges'],
            'colormap': 'bwr',  # 0 = blue, 1 = red
            'width': cellmap.edge_df['edge_width'],
            'alpha': 1.0,
            'head_width': 0.0
        },
        'vert': {
            'visible': True,
            's': 50,
            'color': 'blue',
            'alpha': 1.0
        }
    }

    # Plot
    fig, ax = sheet_view(cellmap, mode="2D", **specs)
    fig.set_size_inches(15, 15)
    ax.set_xticks([])
    ax.set_yticks([])

    return fig, ax


def colour_faces_edges(cellmap, faces, edges):
    """
    Highlights multiple faces and edges on a Tyssue sheet.

    Parameters:
    - cellmap: Tyssue sheet object
    - faces: list of face indices to highlight
    - edges: list of edge indices to highlight
    """
    cellmap.face_df["color_faces"] = 0
    cellmap.edge_df["color_edges"] = 0
    cellmap.edge_df['edge_width'] = 1

    # Highlight given edges
    cellmap.edge_df.loc[edges, "color_edges"] = 1
    cellmap.edge_df.loc[edges, 'edge_width'] = 4

    # Highlight given faces
    cellmap.face_df.loc[faces, "color_faces"] = 1

    specs = {
        'face': {
            'visible': True,
            'color': cellmap.face_df['color_faces'],
            'colormap': 'PuRd',
        },
        'edge': {
            'visible': True,
            'color': cellmap.edge_df['color_edges'],
            'colormap': 'bwr',
            'width': cellmap.edge_df['edge_width'],
        },
        'vert': {
            'visible': False,
        }
    }

    fig, ax = sheet_view(cellmap, **specs)
    fig.set_size_inches((15, 15))


def laser_ablation(cellmap):

    edges_without_opposites_indices, border_edges_indices = cellDivision.edges_list_without_opposites(cellmap)

    chosen_edge = random.choice(edges_without_opposites_indices)

    srce, trgt = cellmap.edge_df.loc[chosen_edge, ["srce", "trgt"]]
    opposite_edge = cellmap.edge_df[(cellmap.edge_df["srce"] == trgt) & (cellmap.edge_df["trgt"] == srce)].index
    opposite_edge = cellmap.edge_df[(cellmap.edge_df["srce"] == trgt) & (cellmap.edge_df["trgt"] == srce)].index.tolist()[0]

    both_edges = [chosen_edge, opposite_edge]

    cellmap_H.edge_df["color_edges"] = 0
    cellmap_H.edge_df['edge_width'] = 1

    cellmap_H.edge_df.loc[both_edges, "color_edges"] = 1
    cellmap_H.edge_df.loc[both_edges, 'edge_width'] = 4

    specs = {
        'face': {
            'visible': False,
        },
        'edge': {
            'visible': True,
            'color': cellmap.edge_df["color_edges"],
            'colormap': 'bwr',
            'width': cellmap.edge_df['edge_width'],
        },
        'vert': {
            'visible': False,
            'color': '#000a4b',
            's': 20,
        }
    }

    fig, ax = sheet_view(cellmap, **specs)
    fig.set_size_inches((15, 15))

    fig, ax = sheet_view(
        cellmap_H,
        mode="2D",
        face={"visible": False},
        edge={"head_width": 0.0, "width": 2},
        vert={"visible": False}
    )

    fig.set_size_inches(20, 20)
    ax.set_xticks([])
    ax.set_yticks([])

    fig, ax

    faces_ids = cellmap.edge_df.loc[both_edges]['face'].tolist()

    srce_vert = cellmap.edge_df.loc[chosen_edge, 'srce']
    trgt_vert = cellmap.edge_df.loc[chosen_edge, 'trgt']
    both_vertices = [srce_vert, trgt_vert]

    initial_coords_srce_vert = cellmap.vert_df.loc[srce_vert, ['x', 'y']]
    initial_coords_trgt_vert = cellmap.vert_df.loc[trgt_vert, ['x', 'y']]

    initial_coords_srce_vert['x', 'y'] = cellmap.vert_df.loc[srce_vert, ['x', 'y']]
    initial_coords_trgt_vert['x', 'y'] = cellmap.vert_df.loc[trgt_vert, ['x', 'y']]

    cellmap.edge_df[(cellmap.edge_df['face'] == faces_ids[0])]

    cellmap.edge_df[(cellmap.edge_df['face'] == faces_ids[1])]

    cellmap.edge_df.loc[cellmap.edge_df['face'] == faces_ids[0], 'face'] = faces_ids[1]

    cellmap.face_df = cellmap.face_df.drop(faces_ids[0])

    cellmap.edge_df = cellmap.edge_df.drop(both_edges)

    cellmap.reset_index()

    energyContributions_model.compute_energy(cellmap)
    [cellmap_H, geom, model_H, history_H, solver1] = vertexModel.solveEuler(cellmap_H, geom, energyContributions_model,
                                                                            50)

    # Store the distances for plotting
    displacement = []
    time_steps = []

    initial_distance = sqrt(
        (initial_coords_srce_vert['x'] - initial_coords_trgt_vert['x']) ** 2 +
        (initial_coords_srce_vert['y'] - initial_coords_trgt_vert['y']) ** 2
    )

    output_folder = "Automated function"
    os.makedirs(output_folder, exist_ok=True)

    file_name = f"FE = {area_elasticity}LE = {length_elasticity}.png"
    file_path = os.path.join(output_folder, file_name)

    displacement.append(0)  # Initial extra distance is zero
    time_steps.append(0)

    for t, cellmap in solver1.history:
        # Extract current coordinates of the source and target vertices
        current_coords_srce_vert = cellmap.vert_df.loc[srce_vert, ['x', 'y']]
        current_coords_trgt_vert = cellmap.vert_df.loc[trgt_vert, ['x', 'y']]

        # Calculate the distance between the source and target vertices
        distance = sqrt(
            (current_coords_srce_vert['x'] - current_coords_trgt_vert['x']) ** 2 +
            (current_coords_srce_vert['y'] - current_coords_trgt_vert['y']) ** 2
        )

        # Calculate the extra distance
        extra_distance = distance - initial_distance

        displacement.append(extra_distance)
        time_steps.append(t)

    # Plotting the extra distances
    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, displacement, marker='o', linestyle='-', color='skyblue',
             label='Extra Distance between Source and Target Vertices')

    plt.title('Displacement Over Time', fontsize=20, fontweight='bold')
    plt.xlabel('Time', fontsize=10, fontweight='bold')
    plt.ylabel('Displacement', fontsize=10, fontweight='bold')
    plt.legend()
    plt.grid(True)
    plt.savefig(file_path)
    plt.show()
    

    return cellmap_H, solver1


def identify_boundary_layers(cellmap, max_layers):
    # Step 1: Add a column 'is_included' to track which faces are included initially
    cellmap.face_df['is_included'] = True  # Initialize all faces as included

    # Keep track of the boundary faces, inside faces, boundary edges, inside edges, inside vertices
    boundary_faces = []
    inside_faces = []
    boundary_edges = []
    inside_edges = []
    inside_vertices = []

    # Keep track of outside vertices
    outside_vertices = set()

    # Keep a working DataFrame for included faces and edges
    working_faces_df = cellmap.face_df.copy()
    working_edges_df = cellmap.edge_df.copy()

    for layer in range(max_layers):
        current_edges_without_opposites = []

        # Step 2: Identify edges without opposites in the current working set
        for edge_idx, edge_row in working_edges_df.iterrows():
            if edge_row['face'] in working_faces_df.index:
                srce, trgt = edge_row['srce'], edge_row['trgt']
                opposite_edges = working_edges_df[
                    (working_edges_df['srce'] == trgt) &
                    (working_edges_df['trgt'] == srce)
                    ]
                if opposite_edges.empty:
                    current_edges_without_opposites.append(edge_idx)

        # Step 3: Exclude faces associated with these edges
        faces_to_exclude = working_edges_df.loc[current_edges_without_opposites, 'face'].unique()

        # Step 4: Exclude other edges associated with the given faces
        other_edges_to_exclude = working_edges_df[working_edges_df['face'].isin(faces_to_exclude)].index.tolist()

        # Add these faces and edges to the boundary list
        boundary_faces.extend(faces_to_exclude)
        boundary_edges.extend(current_edges_without_opposites)
        boundary_edges.extend(other_edges_to_exclude)  # Exclude all edges related to these faces

        # Exclude these faces from the working DataFrame
        working_faces_df = working_faces_df[~working_faces_df.index.isin(faces_to_exclude)]
        working_edges_df = working_edges_df[~working_edges_df['face'].isin(faces_to_exclude)]

        # If no more faces can be excluded, stop the loop
        if working_faces_df.empty:
            break

    # Step 5: Identify inside edges and faces
    inside_faces = working_faces_df.index.tolist()  # The remaining faces
    inside_edges = working_edges_df.index.tolist()  # The remaining edges

    # Step 6: Identify inside vertices by looping over the inside edges
    for edge_id in inside_edges:
        edge = cellmap.edge_df.loc[edge_id]
        srce_vertex = edge['srce']
        trgt_vertex = edge['trgt']

        # Add the vertices to the inside vertices list
        inside_vertices.append(srce_vertex)
        inside_vertices.append(trgt_vertex)

    # Remove duplicates in the inside vertices list
    inside_vertices = list(set(inside_vertices))

    # Step 7: Identify outside vertices (those not in inside vertices)
    for edge_id in boundary_edges:
        edge = cellmap.edge_df.loc[edge_id]
        srce_vertex = edge['srce']
        trgt_vertex = edge['trgt']

        # Add boundary vertices to the outside vertices list if they are not inside vertices
        if srce_vertex not in inside_vertices:
            outside_vertices.add(srce_vertex)
        if trgt_vertex not in inside_vertices:
            outside_vertices.add(trgt_vertex)

    # Convert the set of outside vertices back to a list
    outside_vertices = list(outside_vertices)

    # === NEW: Step 8: Identify outside edges (edges touching outside vertices) ===
    outside_edges = []
    for edge_id, edge in cellmap.edge_df.iterrows():
        if edge['srce'] in outside_vertices or edge['trgt'] in outside_vertices:
            outside_edges.append(edge_id)

    return boundary_edges, boundary_faces, inside_edges, outside_edges, inside_faces, inside_vertices, outside_vertices


def view(cellmap, geom, show_axes=True, xlabel="x", ylabel="y", xlim=None, ylim=None):
    """
    Visualizes the cell map with optional axes and labels.

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


def highlight_vertices(cellmap, geom, chosen_vert_ids, default_color='red', highlight_color='lime', default_size=10, highlight_size=40):
    """
    Visualises the cell map with solid blue edges and highlights selected vertices.

    Parameters:
    - cellmap: Tyssue sheet (cellmap).
    - geom: Geometry object to update cellmap.
    - chosen_vert_ids: List of vertex IDs to highlight.
    - default_color: Default vertex colour (others).
    - highlight_color: Colour for highlighted vertices.
    - default_size: Size for all non-highlighted vertices.
    - highlight_size: Size for highlighted vertices.
    """
    # Update geometry before plotting
    geom.update_all(cellmap)

    # Set default vertex appearance
    cellmap.vert_df['color'] = default_color
    cellmap.vert_df['size'] = default_size

    # Highlight specified vertices
    for v_id in chosen_vert_ids:
        if v_id in cellmap.vert_df.index:
            cellmap.vert_df.loc[v_id, 'color'] = highlight_color
            cellmap.vert_df.loc[v_id, 'size'] = highlight_size

    # Force full opacity for all vertices
    cellmap.vert_df['alpha'] = 1.0

    # Plot with consistent visual style
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
    ax.set_xticks([])
    ax.set_yticks([])

    return fig, ax



def load_simulation_state(folder_path, filename='final_cellmap_state.pkl'):
    file_path = os.path.join(folder_path, filename)

    # Load the cellmap object using pickle
    with open(file_path, 'rb') as f:
        cellmap = pickle.load(f)

    print(f"Simulation state loaded from {file_path}")
    return cellmap


def save_simulation_state(cellmap, folder_path, filename='final_cellmap_state.pkl'):
    # Ensure the folder exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Construct the full file path
    file_path = os.path.join(folder_path, filename)

    # Save the cellmap object using pickle
    with open(file_path, 'wb') as f:
        pickle.dump(cellmap, f)

    print(f"Simulation state saved at {file_path}")




# Register the custom colormap
custom_cmap = ListedColormap(["yellow", "red"])
# Register custom colormap only if it doesn't already exist
if 'yellow_red' not in plt.colormaps():
    custom_cmap = ListedColormap(["yellow", "red"])
    try:
        # Matplotlib 3.5+ preferred API
        plt.colormaps.register(custom_cmap, name='yellow_red')
    except AttributeError:
        # Fallback for older matplotlib versions
        cm.register_cmap(name='yellow_red', cmap=custom_cmap)

def highlight_edge_on_cellmap(cellmap_init, edge_id=None, figsize=(15, 15),
                              vert_col='orange', edge_highlight_col='red',
                              base_edge_col='black', base_vert_col='black', vert_size=240,
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

    from tyssue.draw.plt_draw import draw_edge, draw_vert
    import matplotlib.pyplot as plt

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

    # Draw overlay coloured edges
    draw_edge(
        cellmap,
        coords=coords,
        ax=ax,
        color="yellow",
        #colormap='yellow_red',
        width=cellmap.edge_df["edge_width"],
        head_width=0.0,
        alpha=1.0,
        zorder=2
    )

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
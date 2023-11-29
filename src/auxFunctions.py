import tyssue
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


def length_elasticity_range(cellmap, lower_FRC_length_elasticity, higher_FRC_length_elasticity, lower_ECM_length_elasticity, higher_ECM_length_elasticity):

    cellmap.edge_df['length_elasticity_ECM'] = cellmap.edge_df['length_elasticity']
    cellmap.edge_df['length_elasticity_FRC'] = cellmap.edge_df['length_elasticity']


    for edge in range(len(cellmap.edge_df)):
        newFRCValue = random.uniform(lower_FRC_length_elasticity, higher_FRC_length_elasticity)
        newECMValue = random.uniform(lower_ECM_length_elasticity, higher_ECM_length_elasticity)
        
        
        cellmap.edge_df.loc[edge, 'length_elasticity'] = newFRCValue + newECMValue
        cellmap.edge_df.loc[edge, 'length_elasticity_ECM'] = newECMValue
        cellmap.edge_df.loc[edge, 'length_elasticity_FRC'] = newFRCValue

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

def prefered_area_range(cellmap, a, loc, scale, size):
    newValues = skewnorm.rvs(a, loc, scale, size)  
    positive_values = np.abs(newValues)
    cellmap.face_df['prefered_area'] = np.copy(positive_values)
    
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

    norm = plt.Normalize(vmin=parameter.min(), vmax=parameter.max())
    sm = ScalarMappable(cmap=cm.RdPu, norm=norm)
    sm.set_array([])  # You can set an empty array or a list of values here

    # Add a colorbar using the ScalarMappable
    cbar = plt.colorbar(sm, ax=ax)
    fig.set_size_inches(20,20)

def colour_face(cellmap, face):
    
    cellmap.face_df["color_faces"] = 0
    cellmap.edge_df["color_edges"] = 0
    cellmap.edge_df['edge_width'] = 1


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
    fig.set_size_inches((15, 15))

def colour_edge(cellmap, edge):

    cellmap.edge_df["color_edges"] = 0
    cellmap.edge_df['edge_width'] = 1

    cellmap.edge_df.loc[edge, "color_edges"] = 1
    cellmap.edge_df.loc[edge, 'edge_width'] = 4
    
    specs = {
        'face': {
            'visible': False,
        },
        'edge': {
            'visible': True,
            'color': cellmap.edge_df['color_edges'],
            'colormap':'bwr',
            'width' : cellmap.edge_df['edge_width'],
        },
        'vert': {
            'visible': False,
        }
    }

    fig, ax = sheet_view(cellmap, **specs)
    fig.set_size_inches((15, 15))

def colour_face_edge(cellmap, face, edge):
    
    cellmap.face_df["color_faces"] = 0
    cellmap.edge_df["color_edges"] = 0
    cellmap.edge_df['edge_width'] = 1
    
    cellmap.edge_df.loc[edge, "color_edges"] = 1
    cellmap.edge_df.loc[edge, 'edge_width'] = 4

    cellmap.face_df.loc[face, "color_faces"] = 1
    specs = {
        'face': {
            'visible': True,
            'color': cellmap.face_df['color_faces'],
            'colormap':'PuRd',
        },
        'edge': {
            'visible': True,
            'color': cellmap.edge_df['color_edges'],
            'colormap':'bwr',
            'width' : cellmap.edge_df['edge_width'],
        },
        'vert': {
            'visible': False,
        }
    }

    fig, ax = sheet_view(cellmap, **specs)
    fig.set_size_inches((15, 15))



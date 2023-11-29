import ipyvolume as ipv
import tyssue
import json
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import random
import sys
import tyssue.io.hdf5 as hdf5
from IPython.display import Image
from scipy import optimize
from tyssue.draw import sheet_view
from tyssue.geometry.planar_geometry import PlanarGeometry as geom

##### Own functions
import src.vertexModel as vertexModel
import src.inputMechanicalParameters as inputMechanicalParameters
import src.auxFunctions as auxFunctions
from tyssue.topology import add_vert
import src.cellDivision as cellDivision


from tyssue.topology.base_topology import add_vert, close_face, collapse_edge, remove_face
from tyssue.topology.base_topology import split_vert as base_split_vert

#### Initialize Model
[cellmap_init, geom, energyContributions_model] = vertexModel.initialize()

## Update mechanical parameters
cellmap_init = inputMechanicalParameters.update(cellmap_init)

## Initial stage
energyContributions_model.compute_energy(cellmap_init)

## RUN
[cellmap_H, geom, model_H, history_H] = vertexModel.solveEuler(cellmap_init, geom, energyContributions_model, endTime = 100)


fig, ax = sheet_view(cellmap_H, ['y', 'x'], edge={"color":1, 'colormap': 'Greys'})
#auxFunctions.create_frames(history_H, './results', edge={'color':'black'})
#auxFunctions.exportToMesh(history_H, './results')

import warnings
warnings.filterwarnings("ignore")

### match length_elasticity of opposite edges

for edge in range(len(cellmap_H.edge_df)):

    ### find the opposite edge of each edge in edges_list_without_opposites

    opposite_edge_indices = cellmap_H.edge_df[
        (cellmap_H.edge_df["srce"] == cellmap_H.edge_df.loc[edge, "trgt"]) &
        (cellmap_H.edge_df["trgt"] == cellmap_H.edge_df.loc[edge, "srce"])
        ].index

    if len(opposite_edge_indices):
        # Update the values for the opposite edges
        cellmap_H.edge_df['length_elasticity_ECM'][opposite_edge_indices] = cellmap_H.edge_df['length_elasticity_ECM'][
            edge]
        cellmap_H.edge_df['length_elasticity_FRC'][opposite_edge_indices] = cellmap_H.edge_df['length_elasticity_FRC'][
            edge]
        cellmap_H.edge_df['length_elasticity'][opposite_edge_indices] = cellmap_H.edge_df['length_elasticity'][edge]

### match prefered_length of opposite edges


for edge in range(len(cellmap_H.edge_df)):

    ### find the opposite edge of each edge in edges_list_without_opposites

    opposite_edge_indices = cellmap_H.edge_df[
        (cellmap_H.edge_df["srce"] == cellmap_H.edge_df.loc[edge, "trgt"]) &
        (cellmap_H.edge_df["trgt"] == cellmap_H.edge_df.loc[edge, "srce"])
        ].index

    if len(opposite_edge_indices):
        # Update the values for the opposite edges
        cellmap_H.edge_df['prefered_length'][opposite_edge_indices] = cellmap_H.edge_df['prefered_length'][edge]


def apoptosis_one_FRC(cellmap, condition):
    # Find non-border edges
    edges_without_opposites_indices, border_edges_indices = cellDivision.edges_list_without_opposites(cellmap)

    # Choose edge to remove from non-border edges
    edge_to_remove = cellmap.edge_df[~cellmap.edge_df.index.isin(border_edges_indices)][condition].nlargest(
        1).index.tolist()

    srce, trgt = cellmap.edge_df.loc[edge_to_remove[0], ["srce", "trgt"]]
    opposites = cellmap.edge_df[
        (cellmap.edge_df["srce"] == trgt) & (cellmap.edge_df["trgt"] == srce)]
    opp_edge = opposites.index

    # Remove the edge_to_remove from edge_df
    cellmap.edge_df = cellmap.edge_df[~cellmap.edge_df.index.isin(edge_to_remove)]
    cellmap.edge_df = cellmap.edge_df[~cellmap.edge_df.index.isin(opp_edge)]

    # Reset index and topology (adjust as needed)
    cellmap.reset_index()
    cellmap.reset_topo()

    # Update geometry
    geom.update_all(cellmap)

    energyContributions_model.compute_energy(sheet)
    vertexModel.solveEuler(sheet, geom, energyContributions_model, 40)

    [cellmap_H, geom, model_H, history_H] = vertexModel.solveEuler(cellmap_H, geom, energyContributions_model, endTime = 100)
    auxFunctions.create_frames(history_H, './results', edge={'color': 'black'})
    auxFunctions.exportToMesh(history_H, './results')

    return edge_to_remove

apoptosis_one_FRC(cellmap_H, "length")


import src.auxFunctions as auxFunctions
import numpy as np
from scipy.stats import lognorm


def update(cellmap):
    ## Per sheet:
    cellmap.settings["temperature"] = 10
    #Stochastically detaches vertices from rosettes.
    #Uses two probabilities `p_4` and `p_5p`
    cellmap.settings["p_4"] = 10
    cellmap.settings["p_5p"] = 0.1
    cellmap.settings["threshold_length"] = 2e-2

    ## Per face:

    rangeAreaElasticity = False

    if rangeAreaElasticity:
        lower_area_elasticity = 1
        higher_area_elasticity = 10
        cellmap = auxFunctions.area_elasticity_range(cellmap, lower_area_elasticity, higher_area_elasticity)
    else:
        cellmap.face_df["area_elasticity"] = 1


    #Preferred area is set to mean area
    cellmap.face_df["prefered_area"] = cellmap.face_df["area"].mean()

    ## Per edge, ECM:
    rangeLengthElasticity = False

    if rangeLengthElasticity:
        lower_ECM_le = 12
        higher_ECM_le = 320
        cellmap = auxFunctions.length_elasticity_ECM_range(cellmap, lower_ECM_le, higher_ECM_le)

    else:
        cellmap.edge_df['length_elasticity'] = 100


    ## Per edge, preferred length: 

    rangePreferedLength = False

    if rangePreferedLength:
        lower_prefered_length = 0.05
        higher_prefered_length = 0.1
        cellmap = auxFunctions.prefered_length_range(cellmap, lower_prefered_length, higher_prefered_length)
    else:
        #cellmap.edge_df['prefered_length'] = cellmap.edge_df['length'].mean()
        cellmap.edge_df['prefered_length'] = 0.01


    # Assign FRC tension values if requested
    rangeLineTension = False

    if rangeLineTension:
        lower_FRC_t = 400
        higher_FRC_t = 600
        cellmap = auxFunctions.line_tension_FRC_range(cellmap, lower_FRC_t, higher_FRC_t)

    else:
        cellmap.edge_df["line_tension"] = 380

    ## Per vertex:
    cellmap.vert_df["viscosity"] = 1334
    ## return
    return cellmap
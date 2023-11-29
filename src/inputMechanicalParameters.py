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

    ## Per cell:

    ## Per face:

    #cellmap.face_df["prefered_area"] = cellmap.face_df["area"].mean()*1.1
    #cellmap.face_df["perimeter"] = 1
    #cellmap.face_df["perimeter_elasticity"] = 10
    #cellmap.face_df["prefered_perimeter"] = 3.81

    rangeAreaElasticity = False

    if rangeAreaElasticity:
        lower_area_elasticity = 1
        higher_area_elasticity = 10
        cellmap = auxFunctions.area_elasticity_range(cellmap, lower_area_elasticity, higher_area_elasticity)
    else:
        cellmap.face_df["area_elasticity"] = 0.05


    rangePreferedArea = False

    if rangePreferedArea:
        size = len(cellmap.face_df)
        loc = 6
        scale = 1
        a = 100
        cellmap = auxFunctions.prefered_area_range(cellmap, a, loc, scale, size)
    else:
        cellmap.face_df["prefered_area"] = 10


    ## Per edge:

    

    rangeLengthElasticity = True

    if rangeLengthElasticity:
        lower_FRC_length_elasticity = 80
        higher_FRC_length_elasticity = 1200
        lower_ECM_length_elasticity = 40
        higher_ECM_length_elasticity = 60
        cellmap = auxFunctions.length_elasticity_range(cellmap, lower_FRC_length_elasticity, higher_FRC_length_elasticity, lower_ECM_length_elasticity, higher_ECM_length_elasticity)
    else:
        cellmap.edge_df["length_elasticity"] = 100

    cellmap = auxFunctions.set_opposite_edges_mechanics(cellmap)


    rangePreferedLength = True

    if rangePreferedLength:
        lower_prefered_length = 0.05
        higher_prefered_length = 1
        cellmap = auxFunctions.prefered_length_range(cellmap, lower_prefered_length, higher_prefered_length)
    else:
        cellmap.edge_df["prefered_length"] = 0.05

    ## Per vertex:
    cellmap.vert_df["viscosity"] = 1

    ## return
    return cellmap
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
        cellmap.face_df["area_elasticity"] = 10


    rangePreferedArea = False
    splitPreferedArea = False

    # Set preferred face areas (FRC areas)
    if rangePreferedArea:
        # Apply a range for preferred area to all FRCs
        lower_prefered_area = 8.13
        higher_prefered_area = 11.87
        cellmap = auxFunctions.prefered_area_range(cellmap, lower_prefered_area, higher_prefered_area)

    elif splitPreferedArea:
        # Split preferred area into two distinct values
        area_x = 10  # Preferred area for the first half of FRCs
        area_y = 40  # Preferred area for the second half of FRCs
        cellmap = auxFunctions.random_island_prefered_area(cellmap, area_x, area_y)
    else:
        # Default preferred area if neither range nor split is specified
        cellmap.face_df["prefered_area"] = cellmap.face_df["area"].mean()
        #cellmap.face_df["prefered_area"] = 50


        ## Per edge:
    rangeLengthElasticity = False

    if rangeLengthElasticity:
        lower_FRC_le = 11
        higher_FRC_le = 18
        lower_ECM_le = 5
        higher_ECM_le = 145
        cellmap = auxFunctions.length_elasticity_range(cellmap, lower_FRC_le, higher_FRC_le,
                                                       lower_ECM_le, higher_ECM_le)
    else:
        cellmap.edge_df["length_elasticity"] = 90

    #cellmap = auxFunctions.set_opposite_edges_mechanics(cellmap)

    rangePreferedLength = False

    if rangePreferedLength:
        lower_prefered_length = 0.05
        higher_prefered_length = 0.1
        cellmap = auxFunctions.prefered_length_range(cellmap, lower_prefered_length, higher_prefered_length)
    else:
        cellmap.edge_df['prefered_length'] = 1
        #cellmap.edge_df['prefered_length'] = cellmap.edge_df['length'].mean()

    ## Per vertex:
    cellmap.vert_df["viscosity"] = 1268
    ## return
    return cellmap
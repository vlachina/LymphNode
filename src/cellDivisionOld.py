import src.auxFunctions as auxFunctions
import src.vertexModel as vertexModel 
import logging 
import pandas as pd 
from tyssue import Sheet, PlanarGeometry  as geom 
from tyssue.solvers.viscous import EulerSolver 
 
 
 
def ln_divide_FRC(eptm, edge): 
 
    ### function adapted from add_vert in tyssue (https://tyssue.readthedocs.io/en/latest/source/tyssue.topology.html#tyssue.topology.base_topology.add_vert) 
 
    srce, trgt = eptm.edge_df.loc[edge, ["srce", "trgt"]] 
 
    logger = logging.getLogger(name=__name__) 
    logger.debug(f"adding vertex between {srce} and {trgt}") 
     
   ### find the opposite half-edge 
 
    opposites = eptm.edge_df[ 
    (eptm.edge_df["srce"] == trgt) & (eptm.edge_df["trgt"] == srce)] 
    opp_edges = opposites.index 
 
    ### find any parallel edges 
 
    parallels = eptm.edge_df[ 
    (eptm.edge_df["srce"] == srce) & (eptm.edge_df["trgt"] == trgt)]

    ### creates a new row in vert_df that has a specific srce value
    new_vert = eptm.vert_df.loc[srce:srce]
     
   ### joins the new new_vert row with vert_df 
 
    eptm.vert_df = pd.concat([eptm.vert_df, new_vert], ignore_index=True) 
     
   ### assigns an index to the new_vert at the end of dataframe 
 
    new_vert = eptm.vert_df.index[-1] 
     
   ### updates its coordinates with mean of original srce and trgt vertices 
 
    eptm.vert_df.loc[new_vert, eptm.coords] = eptm.vert_df.loc[ 
    [srce, trgt], eptm.coords 
    ].mean(numeric_only=True) 
 
    eptm.edge_df.loc[parallels.index, "trgt"] = new_vert 
     
   ### add new edges 
 
    eptm.edge_df = pd.concat([eptm.edge_df, parallels], ignore_index=True) 
    new_edges = eptm.edge_df.index[-parallels.index.size :] 
     
   ### update src and trg of new_edges 
 
    eptm.edge_df.loc[new_edges, "srce"] = new_vert 
    eptm.edge_df.loc[new_edges, "trgt"] = trgt 
     
   ### add new opposite edges 
 
    new_opp_edges = [] 
 
    if len(opposites.index): 
        eptm.edge_df.loc[opposites.index, "srce"] = new_vert 
        eptm.edge_df = pd.concat([eptm.edge_df, opposites], ignore_index=True) 
        new_opp_edges = eptm.edge_df.index[-opposites.index.size :] 
         
   ### assign trgt and srce to new_opposite_edges 
 
        eptm.edge_df.loc[new_opp_edges, "trgt"] = new_vert 
        eptm.edge_df.loc[new_opp_edges, "srce"] = trgt 
 
 
    new_length_elasticity_FRC = eptm.edge_df.loc[edge, 'length_elasticity_FRC'] / 1.5 
 
    eptm.edge_df.loc[edge, 'length_elasticity_FRC'] = new_length_elasticity_FRC 
    eptm.edge_df.loc[new_edges, "length_elasticity_FRC"] = new_length_elasticity_FRC 
    eptm.edge_df.loc[new_opp_edges, "length_elasticity_FRC"] = new_length_elasticity_FRC 
    eptm.edge_df.loc[opp_edges, 'length_elasticity_FRC'] = new_length_elasticity_FRC 
 
 
    new_prefered_length_FRC = eptm.edge_df.loc[edge, 'prefered_length'] / 2 
 
    eptm.edge_df.loc[edge, 'prefered_length'] = new_prefered_length_FRC 
    eptm.edge_df.loc[new_edges, "prefered_length"] = new_prefered_length_FRC 
    eptm.edge_df.loc[new_opp_edges, "prefered_length"] = new_prefered_length_FRC 
    eptm.edge_df.loc[opp_edges, 'prefered_length'] = new_prefered_length_FRC 
 
 
    # ## Sheet special case 
    if len(new_edges) == 1: 
        new_edges = new_edges[0] 
    if len(new_opp_edges) == 1: 
        new_opp_edges = new_opp_edges[0] 
    elif len(new_opp_edges) == 0: 
        new_opp_edges = None 
    return new_vert, new_edges, new_opp_edges, opp_edges 
 
 
 
### make a list of edges without their opposite edges or border edges 
 
 
def edges_list_without_opposites(cellmap):
    # Create a set to keep track of edge pairs that have been processed
    edges_without_opposites_indices = []
    border_edges_indices = []
    processed_edge_pairs = set() 
 
    for idx, edge in cellmap.edge_df.iterrows(): 
        srce, trgt = edge['srce'], edge['trgt'] 
        reversed_edge = (cellmap.edge_df['srce'] == trgt) & (cellmap.edge_df['trgt'] == srce) 
        edge_pair_id = tuple(sorted((srce, trgt))) 
 
    # Check if the edge pair has already been processed and it's an edge without an opposite edge 
        if edge_pair_id not in processed_edge_pairs and reversed_edge.any(): 
        # Mark the edge pair as processed 
            processed_edge_pairs.add(edge_pair_id) 
        # Add the edge to the filtered list 
            edges_without_opposites_indices.append(idx)


        if not reversed_edge.any():
            border_edges_indices.append(idx)
 
    return edges_without_opposites_indices, border_edges_indices 
 
 
def ln_condition_for_division_highest(sheet, condition, percent):

    edges_without_opposites_indices, border_edges_indices = edges_list_without_opposites(sheet)
 
    filtered_edges = sheet.edge_df.loc[edges_without_opposites_indices]
    num_to_select = int(percent / 100 * len(filtered_edges)) 
    chosen_edges = filtered_edges.sort_values(by=condition, ascending=False).head(num_to_select)

    chosen_edges_indices = chosen_edges.index.tolist()
    return chosen_edges_indices 
 
 
def ln_condition_for_division_lowest(sheet, condition, percent):

    edges_without_opposites_indices, border_edges_indices = edges_list_without_opposites(sheet)
 
    filtered_edges = sheet.edge_df.loc[edges_without_opposites_indices] 
    num_to_select = int(percent / 100 * len(filtered_edges)) 
    chosen_edges = filtered_edges.sort_values(by=condition, ascending=True).head(num_to_select)

    chosen_edges_indices = chosen_edges.index.tolist()
    return chosen_edges_indices
 
def ln_divide_FRCs(sheet, condition, percent, energyContributions_model):

    new_edges_indices = []
    new_opp_edges_indices = []
    new_vertices_indices = []

    chosen_edges_indices = ln_condition_for_division_highest(sheet, condition, percent)
  
    for edge in chosen_edges_indices:

        print("Chosen edge:", edge)
         
        new_vert, new_edges, new_opp_edges, opp_edges = ln_divide_FRC(sheet, edge)
        new_edges_indices.append(new_edges)
        new_opp_edges_indices.append(new_opp_edges)
        new_vertices_indices.append(new_vert)
 
         
        print("New vert:", new_vert)
        print("New Edge:", new_edges)
        print("New opposite edge:", new_opp_edges)
                 
        geom.update_all(sheet)
         
        energyContributions_model.compute_energy(sheet)
        vertexModel.solveEuler(sheet, geom, energyContributions_model, 20)
 
 
    return sheet, new_edges_indices, new_opp_edges_indices, new_vertices_indices 
 
 

import random


def find_vertices_of_face(sheet, face):
    # Get the edges associated with the given face
    face_edges = sheet.edge_df[sheet.edge_df['face'] == face].index.tolist()

    # Extract unique vertices from the edges
    vertices = set()
    for edge_id in face_edges:
        edge = sheet.edge_df.loc[edge_id]
        vertices.add(edge['srce'])
        vertices.add(edge['trgt'])

    return list(vertices)


def find_edges_a_and_b(sheet, mother):

    ### find edge_a which is the new edge created during ln_divide_FRC
    edges_a = cellmap_H.edge_df.loc[new_edges_indices].index.tolist()

    ### find edge_a for each mother cell
    print("Mother cell: ", mother)
    mother_cell_edges = cellmap_H.edge_df[cellmap_H.edge_df['face'] == mother].index.tolist()

    edge_a = None
    edge_b = None

    for edge in mother_cell_edges:
        if edge in edges_a:
            edge_a = edge
            print("Chosen edge_a:", edge_a)
    if edge_a is None:
        print("No edge_a found in the given mother cell.")

    ### find edge_b for each mother cell

    #exclude edge_a from choices for edge_b
    remaining_edges = [edge for edge in mother_cell_edges if edge not in edges_a]

    #set up a finite amount of attempts to avoid an infinite loop
    attempts = 0

    #check remaining edges aren't empty and the loop has been attempted less than 20 times
    while remaining_edges and attempts < 20:

        # Choose edge_b randomly from the remaining edges
        potential_edge_b = random.choice(remaining_edges)
        print(potential_edge_b)
        #check that edge_b isn't connected to edge_a
        if (cellmap_H.edge_df.loc[potential_edge_b]["srce"] != cellmap_H.edge_df.loc[edge_a]["trgt"]) & (
                cellmap_H.edge_df.loc[potential_edge_b]["trgt"] != cellmap_H.edge_df.loc[edge_a]["srce"]):
            edge_b = potential_edge_b
            print("Chosen edge_b:", edge_b)
            break

        #remove tried edge_b from remaining edges to avoid it being picked again
        remaining_edges.remove(potential_edge_b)
        attempts += 1

    if edge_b is None:
        print("No edge_b found in the given mother cell.")

    return edge_a, edge_b

#function made to connect the protrusion to a random place on edge_b
def ln_divide_FRC_random_angle(sheet, edge):
    ### function adapted from add_vert in tyssue (https://tyssue.readthedocs.io/en/latest/source/tyssue.topology.html#tyssue.topology.base_topology.add_vert)
    """Adds a vertex in the middle of the edge,

    which is split as is its opposite(s)

    Parameters
    ----------
    sheet : a :class:`Epithelium` instance
    edge : int
    the index of one of the half-edges to split

    Returns
    -------
    new_vert : int
    the index to the new vertex
    new_edges : int or list of ints
    index to the new edge(s). For a sheet, returns
    a single index, for a 3D epithelium, returns
    the list of all the new parallel edges
    new_opp_edges : int or list of ints
    index to the new opposite edge(s). For a sheet, returns
    a single index, for a 3D epithelium, returns
    the list of all the new parallel edges


    In the simple case whith two half-edge, returns
    indices to the new edges, with the following convention:

    s    e    t
      ------>
    * <------ *
    oe

    s    e       ne   t
      ------   ----->
    * <----- * ------ *
        oe   nv   noe

    where "e" is the passed edge as argument, "s" its source "t" its
    target and "oe" its opposite. The returned edges are the ones
    between the new vertex and the input edge's original target.
    """

    srce, trgt = sheet.edge_df.loc[edge, ["srce", "trgt"]]

    logger = logging.getLogger(name=__name__)
    logger.debug(f"adding vertex between {srce} and {trgt}")
    opposites = sheet.edge_df[
        (sheet.edge_df["srce"] == trgt) & (sheet.edge_df["trgt"] == srce)
        ]
    parallels = sheet.edge_df[
        (sheet.edge_df["srce"] == srce) & (sheet.edge_df["trgt"] == trgt)
        ]

    new_vert = sheet.vert_df.loc[srce:srce]
    sheet.vert_df = pd.concat([sheet.vert_df, new_vert], ignore_index=True)
    new_vert = sheet.vert_df.index[-1]

    random_weight = random.uniform(0, 1)

    # Calculate the new coordinates as a weighted average of the source and target coordinates
    sheet.vert_df.loc[new_vert, sheet.coords] = (
            random_weight * sheet.vert_df.loc[srce, sheet.coords]
            + (1 - random_weight) * sheet.vert_df.loc[trgt, sheet.coords])

    sheet.edge_df.loc[parallels.index, "trgt"] = new_vert
    sheet.edge_df = pd.concat([sheet.edge_df, parallels], ignore_index=True)
    new_edges = sheet.edge_df.index[-parallels.index.size:]
    sheet.edge_df.loc[new_edges, "srce"] = new_vert
    sheet.edge_df.loc[new_edges, "trgt"] = trgt

    sheet.edge_df.loc[opposites.index, "srce"] = new_vert
    sheet.edge_df = pd.concat([sheet.edge_df, opposites], ignore_index=True)
    new_opp_edges = sheet.edge_df.index[-opposites.index.size:]
    sheet.edge_df.loc[new_opp_edges, "trgt"] = new_vert
    sheet.edge_df.loc[new_opp_edges, "srce"] = trgt

    # ## Sheet special case
    if len(new_edges) == 1:
        new_edges = new_edges[0]
    if len(new_opp_edges) == 1:
        new_opp_edges = new_opp_edges[0]
    elif len(new_opp_edges) == 0:
        new_opp_edges = None
    return new_vert, new_edges, new_opp_edges


def ln_face_division(sheet, mother, vert_a, vert_b):

    ### function adapted from face_division in tyssue

    """
    Divides the face associated with edges
    indexed by `edge_a` and `edge_b`, splitting it
    in the middle of those edes.
    """
    # mother = sheet.edge_df.loc[edge_a, 'face']

    face_cols = sheet.face_df.loc[mother:mother]

    sheet.face_df = pd.concat([sheet.face_df, face_cols], ignore_index=True)
    sheet.face_df.index.name = "face"
    daughter = int(sheet.face_df.index[-1])

    edge_cols = sheet.edge_df[sheet.edge_df["face"] == mother].iloc[0:1]

    sheet.edge_df = pd.concat([sheet.edge_df, edge_cols, edge_cols], ignore_index=True)
    new_edge_m = sheet.edge_df.index[-2]
    sheet.edge_df.loc[new_edge_m, "srce"] = vert_b
    sheet.edge_df.loc[new_edge_m, "trgt"] = vert_a
    new_edge_d = sheet.edge_df.index[-1]
    sheet.edge_df.loc[new_edge_d, "srce"] = vert_a
    sheet.edge_df.loc[new_edge_d, "trgt"] = vert_b

    ### Discover daughter edges
    m_data = sheet.edge_df[sheet.edge_df["face"] == mother]
    daughter_edges = [new_edge_d]
    srce, trgt = vert_a, vert_b
    srces, trgts = m_data[["srce", "trgt"]].values.T
    spins = 0

    while trgt != vert_a:
        srce, trgt = trgt, trgts[srces == trgt][0]

        daughter_edges.append(
            m_data[(m_data["srce"] == srce) & (m_data["trgt"] == trgt)].index[0]
        )
        spins += 1
        if spins > m_data.shape[0]:
            raise ValueError(f"The face {mother} has an invalid topology, \n")
    sheet.edge_df.loc[daughter_edges, "face"] = daughter
    sheet.edge_df.index.name = "edge"
    sheet.reset_topo()

    return daughter, daughter_edges, new_edge_d


def ln_form_protrusion(sheet, mother, geom, angle=None):
    """Causes a cell to divide

    Parameters
    ----------

    sheet : a 'Sheet' instance
    mother : face index of target dividing cell
    geom : a 2D geometry
    angle : division angle for newly formed edge

    Returns
    -------
    daughter: face index of new cell

    Notes
    -----
    - Function checks for perodic boundaries if there are, it checks if dividing cell
      rests on an edge of the periodic boundaries if so, it displaces the boundaries
      by a half a period and moves the target cell in the bulk of the tissue. It then
      performs cell division normally and reverts the periodic boundaries
      to the original configuration
    """

    if sheet.settings.get("boundaries") is not None:
        mother_on_periodic_boundary = False
        if (
                sheet.face_df.loc[mother]["at_x_boundary"]
                or sheet.face_df.loc[mother]["at_y_boundary"]
        ):
            mother_on_periodic_boundary = True
            saved_boundary = sheet.specs["settings"]["boundaries"].copy()
            for u, boundary in sheet.settings["boundaries"].items():
                if sheet.face_df.loc[mother][f"at_{u}_boundary"]:
                    period = boundary[1] - boundary[0]
                    sheet.specs["settings"]["boundaries"][u] = [
                        boundary[0] + period / 2.0,
                        boundary[1] + period / 2.0,
                    ]
            geom.update_all(sheet)

    if not sheet.face_df.loc[mother, "is_alive"]:
        logger.warning("Cell %s is not alive and cannot devide", mother)
        return
    edge_a, edge_b = find_edges_a_and_b(sheet, mother)
    if edge_a is None or edge_b is None:
        print("Skipping cell division due to missing edges.")
        return None, None, None

    ### find vertices associated with mother cell
    vertices_of_mother = find_vertices_of_face(sheet, mother)

    # Find the vertex in new_vertices_indices associated with the mother cell
    vert_a = None
    for vert in new_vertices_indices:
        if vert in vertices_of_mother:
            vert_a = vert
            break
    if vert_a is None:
        print("No vert_a found.")

    vert_b, *_ = ln_divide_FRC_random_angle(sheet, edge_b)
    sheet.vert_df.index.name = "vert"
    daughter, daughter_edges, new_edge_d = face_division(sheet, mother, vert_a, vert_b)

    if sheet.settings.get("boundaries") is not None and mother_on_periodic_boundary:
        sheet.specs["settings"]["boundaries"] = saved_boundary
        geom.update_all(sheet)
    geom.update_all(sheet)

    energyContributions_model.compute_energy(sheet)
    vertexModel.solveEuler(sheet, geom, energyContributions_model, 20)

    return daughter, daughter_edges, new_edge_d






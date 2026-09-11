import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(name=__name__)

def collapse_single_edge_expansion(cellmap, geom, energyContributions_model, edge_id):
    """
    Collapses a single specified edge by merging its two vertices into one.
    The new vertex is placed at the midpoint of the original edge.
    
    Parameters:
    -----------
    cellmap : object
        The cellmap containing vertex and edge DataFrames
    geom : object
        Geometry handler
    energyContributions_model : object
        Energy model for the system
    edge_id : int
        ID of the edge to collapse (must be an inside edge, not boundary)
    
    Returns:
    --------
    cellmap : object
        Updated cellmap after edge collapse
    """
    
    logger.info(f"Collapsing edge: {edge_id}")
    
    # Verifying edge exists
    if edge_id not in cellmap.edge_df.index:
        raise ValueError(f"Edge {edge_id} not found in cellmap.edge_df")
    
    # Finding the two vertices belonging to the edge
    edge_row = cellmap.edge_df.loc[edge_id]
    v1 = edge_row['srce']
    v2 = edge_row['trgt']
    
    logger.info(f"Merging vertices {v1} and {v2}")
    
    # Verifying both vertices exist
    if v1 not in cellmap.vert_df.index or v2 not in cellmap.vert_df.index:
        raise ValueError(f"Vertex {v1} or {v2} not found in cellmap.vert_df")
    
    # Calculating midpoint coordinates along the chosen edge 
    x1 = cellmap.vert_df.loc[v1, 'x']
    y1 = cellmap.vert_df.loc[v1, 'y']
    x2 = cellmap.vert_df.loc[v2, 'x']
    y2 = cellmap.vert_df.loc[v2, 'y']
    
    midpoint_x = (x1 + x2) / 2
    midpoint_y = (y1 + y2) / 2
    
    # Creating a new vertex at the midpoint
    new_vertex_id = max(cellmap.vert_df.index) + 1 if not cellmap.vert_df.empty else 0
    new_vertex_data = cellmap.vert_df.loc[v1].copy()
    new_vertex_data['x'] = midpoint_x
    new_vertex_data['y'] = midpoint_y
    new_vertex_data['new_vert_id'] = np.nan   # ensures new_vert_ids assigns a novel one later on, not v1's old ID
    cellmap.vert_df.loc[new_vertex_id] = new_vertex_data
    
    # Rewiring edges by vertex proximity
    cellmap.edge_df.loc[cellmap.edge_df["srce"] == v1, "srce"] = new_vertex_id
    cellmap.edge_df.loc[cellmap.edge_df["srce"] == v2, "srce"] = new_vertex_id
    cellmap.edge_df.loc[cellmap.edge_df["trgt"] == v1, "trgt"] = new_vertex_id
    cellmap.edge_df.loc[cellmap.edge_df["trgt"] == v2, "trgt"] = new_vertex_id
    
    # Deleting the collapsed edge and its parallel edge
    parallel_edges = cellmap.edge_df[
        ((cellmap.edge_df["srce"] == v1) & (cellmap.edge_df["trgt"] == v2)) |
        ((cellmap.edge_df["srce"] == v2) & (cellmap.edge_df["trgt"] == v1))
    ].index.tolist()
    
    edges_to_delete = [edge_id] + parallel_edges
    cellmap.edge_df.drop(edges_to_delete, inplace=True, errors='ignore')
    
    # Deleting old vertices
    cellmap.vert_df.drop([v1, v2], inplace=True, errors='ignore')
    
    # Removing self-loops (safety mechanism, shouldn't exist)
    cellmap.edge_df = cellmap.edge_df[cellmap.edge_df["srce"] != cellmap.edge_df["trgt"]]
    
    # Removing duplicate edges (if many consecutive collapses form any)
    cellmap.edge_df = cellmap.edge_df.drop_duplicates(subset=['srce', 'trgt'])
    
    # Resetting indices and updating geometry
    cellmap.reset_index()
    
    # Updating active vertices
    if hasattr(cellmap, 'active_verts'):
        cellmap.active_verts = list(cellmap.vert_df.index)
    
    geom.update_all(cellmap)
    
    # Recomputing energy and relaxing the model post edge contraction
    energyContributions_model.compute_energy(cellmap)
    [cellmap, geom, model_H, history_H, solver] = vertexModel2.solveEuler(
        cellmap, geom, energyContributions_model, endTime=40
    )
    
    return cellmap

def split_vertex_expansion(cellmap, chosen_vertex, geom, energyContributions_model, distance, retry_attempts=3):
    """
    Safely divide a vertex by creating a nearby new vertex and rewiring edges.
    On failure, rollback and retry up to `retry_attempts` times.

    Returns:
        (cellmap, chosen_vertex, new_vert_index, new_edge_index, opposite_edge_index)
        or (cellmap, chosen_vertex, None, None, None) if all attempts fail.
    """

    original_vert_df = cellmap.vert_df.copy()
    original_edge_df = cellmap.edge_df.copy()

    new_vert_index = None
    new_edge_index = None
    opposite_edge_index = None

    attempt = 0
    while attempt < retry_attempts:
        try:
            # Creating temorary cellmap copies to allow rollback
            temp_vert_df = cellmap.vert_df.copy()
            temp_edge_df = cellmap.edge_df.copy()

            # Finding connected edges to the chosen vertex
            connected_edges = temp_edge_df[
                (temp_edge_df['srce'] == chosen_vertex) |
                (temp_edge_df['trgt'] == chosen_vertex)
            ].copy()

            # Creating new vertex inside a specified radius from chosen vertex
            new_vert_data = temp_vert_df.loc[chosen_vertex].copy()
            angle = np.random.uniform(0, 2*np.pi)
            dx = distance * np.cos(angle)
            dy = distance * np.sin(angle)
            new_vert_data[cellmap.coords] = temp_vert_df.loc[chosen_vertex, cellmap.coords] + [dx, dy]

            new_vert_index = int(temp_vert_df.index.max()) + 1
            temp_vert_df.loc[new_vert_index] = new_vert_data
    
        
            # Creating a new pair of edges
            # Template: first connected edge (copies its all mechanical properties)

            source_edge = connected_edges.iloc[0]

            template = source_edge.copy()

            new_edge_index = int(temp_edge_df.index.max()) + 1
            opposite_edge_index = new_edge_index + 1

            new_edge = template.copy()
            new_edge["srce"], new_edge["trgt"] = chosen_vertex, new_vert_index
            new_edge["face"] = np.nan

            opposite_edge = template.copy()
            opposite_edge["srce"], opposite_edge["trgt"] = new_vert_index, chosen_vertex
            opposite_edge["face"] = np.nan

            temp_edge_df.loc[new_edge_index] = new_edge
            temp_edge_df.loc[opposite_edge_index] = opposite_edge

            # Reassigning original edges to closer vertex
            chosen_xy = temp_vert_df.loc[chosen_vertex, cellmap.coords].values.astype(float)
            new_xy    = temp_vert_df.loc[new_vert_index, cellmap.coords].values.astype(float)

            reassigned_edges_to_new_vertex = []
            for e_idx, e in connected_edges.iterrows():
                if e_idx in (new_edge_index, opposite_edge_index):
                    continue
                other = e['trgt'] if e['srce'] == chosen_vertex else e['srce']
                other_xy = temp_vert_df.loc[other, cellmap.coords].values.astype(float)

                d_chosen = np.linalg.norm(other_xy - chosen_xy)
                d_new    = np.linalg.norm(other_xy - new_xy)

                if d_new < d_chosen:
                    reassigned_edges_to_new_vertex.append(e_idx)
                    if e['srce'] == chosen_vertex:
                        temp_edge_df.loc[e_idx, 'srce'] = new_vert_index
                    else:
                        temp_edge_df.loc[e_idx, 'trgt'] = new_vert_index

            # Identifying "open" faces (whose edge chains don't close)
            open_faces = []
            for face_id, group in temp_edge_df.groupby('face'):
                verts = list(group[['srce', 'trgt']].itertuples(index=False, name=None))
                if not verts:
                    continue

                # try to follow the loop
                chain = [verts[0][0], verts[0][1]]
                used = {0}
                while True:
                    extended = False
                    for i, (s, t) in enumerate(verts):
                        if i in used:
                            continue
                        if chain[-1] == s:
                            chain.append(t)
                            used.add(i)
                            extended = True
                            break
                        elif chain[-1] == t:
                            chain.append(s)
                            used.add(i)
                            extended = True
                            break
                    if not extended:
                        break

                # if loop doesn't close, this face is "open"
                if chain[0] != chain[-1]:
                    open_faces.append(face_id)

            # filtering to only faces touching the new vertex
            open_faces_touching_new = []
            for f in open_faces:
                verts_f = temp_edge_df[temp_edge_df['face'] == f][['srce', 'trgt']].values.ravel()
                if new_vert_index in verts_f or chosen_vertex in verts_f:
                    open_faces_touching_new.append(f)

            print("Open faces touching new vertex:", open_faces_touching_new)

            if len(open_faces_touching_new) != 2:
                raise ValueError(
                    f"Expected 2 open faces, found {len(open_faces_touching_new)}: {open_faces_touching_new}"
                )

            # Finding which new edge closes which open face (needs correct directionality)
            for f in open_faces_touching_new:
                f_edges = temp_edge_df[temp_edge_df['face'] == f]

                # collecting src/trgt vertices
                srces = list(f_edges['srce'].astype(int))
                trgts = list(f_edges['trgt'].astype(int))

                # imbalance: start and end
                start_candidates = [v for v in srces if v not in trgts]
                end_candidates   = [v for v in trgts if v not in srces]

                if len(start_candidates) == 1 and len(end_candidates) == 1:
                    start = start_candidates[0]
                    end   = end_candidates[0]
                    needed_edge = (end, start)  # must go end→start to close loop

                    new_pair      = (int(temp_edge_df.loc[new_edge_index, 'srce']),
                                     int(temp_edge_df.loc[new_edge_index, 'trgt']))
                    opposite_pair = (int(temp_edge_df.loc[opposite_edge_index, 'srce']),
                                     int(temp_edge_df.loc[opposite_edge_index, 'trgt']))

                    if new_pair == needed_edge:
                        temp_edge_df.loc[new_edge_index, 'face'] = f
                    elif opposite_pair == needed_edge:
                        temp_edge_df.loc[opposite_edge_index, 'face'] = f
                    else:
                        print(f"Face {f}: expected {needed_edge}, "
                              f"but new={new_pair}, opp={opposite_pair}")
                        
            # Verifying both faces are closed in directed edge cycles; if not, swap once and recheck

            def _face_closes(df, face_id):
                sub = df[df['face'] == face_id][['srce', 'trgt']]
                # in==out at every vertex
                outc = sub['srce'].value_counts()
                inc  = sub['trgt'].value_counts()
                verts = set(outc.index) | set(inc.index)
                for v in verts:
                    if outc.get(v, 0) != inc.get(v, 0):
                        return False
                # follow edges as a walk using srce->trgt
                start = int(sub.iloc[0]['srce'])
                cur = start
                used = set()
                for _ in range(len(sub)):
                    nxt = sub[~sub.index.isin(used) & (sub['srce'] == cur)]
                    if nxt.empty:
                        return False
                    eidx = nxt.index[0]
                    used.add(eidx)
                    cur = int(sub.loc[eidx, 'trgt'])
                return cur == start and len(used) == len(sub)

            face_new = int(temp_edge_df.loc[new_edge_index, 'face'])
            face_opp = int(temp_edge_df.loc[opposite_edge_index, 'face'])

            ok_new = _face_closes(temp_edge_df, face_new)
            ok_opp = _face_closes(temp_edge_df, face_opp)

            if not (ok_new and ok_opp):
                # try swapping faces between the two new edges once
                temp_edge_df.loc[new_edge_index, 'face'], temp_edge_df.loc[opposite_edge_index, 'face'] = face_opp, face_new
                face_new, face_opp = face_opp, face_new
                ok_new = _face_closes(temp_edge_df, face_new)
                ok_opp = _face_closes(temp_edge_df, face_opp)

            if not (ok_new and ok_opp):
                raise ValueError(
                    f"Face closure failed after assignment: "
                    f"new→face {face_new} ok={ok_new}, opp→face {face_opp} ok={ok_opp}"
            )

            # Committing temp results
            cellmap.vert_df = temp_vert_df
            cellmap.edge_df = temp_edge_df

            # Updating geometry
            geom.update_all(cellmap)
            cellmap.reset_topo()
            cellmap.reset_index()

            # Relaxing the model
            energyContributions_model.compute_energy(cellmap)
            [cellmap, geom, model_H, history_H, solver] = vertexModel2.solveEuler(
                cellmap, geom, energyContributions_model, endTime=40
            )

            print(f" Successfully divided vertex {chosen_vertex} → new vertex {new_vert_index}")
            return cellmap, chosen_vertex, new_vert_index, new_edge_index, opposite_edge_index

        except Exception as e:
            print(f" split_vertex attempt {attempt+1}/{retry_attempts} failed: {e}")
            # rollback original state
            cellmap.vert_df = original_vert_df.copy()
            cellmap.edge_df = original_edge_df.copy()
            geom.update_all(cellmap)
            cellmap.reset_topo()
            cellmap.reset_index()
            attempt += 1

    print("Failed to divide after multiple attempts.")
    return cellmap, chosen_vertex, None, None, None

def identify_vertices_to_divide(cellmap, inside_vertices, edge_sum_threshold):
    """
    Identify vertices that need division based on edge length sum.
    """
    vertices_to_divide = []
    for v in inside_vertices:
        connected_edges = cellmap.edge_df[
            (cellmap.edge_df["srce"] == v) |
            (cellmap.edge_df["trgt"] == v)
        ].copy()

        # Remove duplicate edges (A→B and B→A)
        connected_edges['pair'] = connected_edges.apply(
            lambda row: tuple(sorted((row['srce'], row['trgt']))), axis=1
        )
        connected_edges = connected_edges.drop_duplicates(subset='pair')

        edge_sum = connected_edges['length'].sum()

        if edge_sum > edge_sum_threshold:
            vertices_to_divide.append(v)
    
    return vertices_to_divide

def perform_divisions(cellmap, vertices_to_divide, next_edge_uid, geom, energyContributions_model, division_distance):
    """
    Perform vertex divisions and return:
      - chosen_vertices: list of vertices that were divided
      - new_vertices: list of newly created vertices
    """
    chosen_vertices = []
    new_vertices = []

    for v_id in vertices_to_divide:
        try:
            cellmap, chosen_vertex, new_vert_index, new_edge_index, opposite_edge_index, next_edge_uid = split_vertex(
                cellmap, v_id, next_edge_uid, geom, energyContributions_model, distance=division_distance
            )

            # Only record if a new vertex was actually created
            if new_vert_index is not None:
                chosen_vertices.append(chosen_vertex)
                new_vertices.append(new_vert_index)

            # Reset after division attempt (safe either way)
            cellmap.reset_index()
            cellmap.reset_topo()
            geom.update_all(cellmap)

        except Exception as e:
            print(f"Division failed for vertex {v_id}: {e}")

    return cellmap, chosen_vertices, new_vertices, next_edge_uid

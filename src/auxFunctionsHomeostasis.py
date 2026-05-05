import numpy as np
import pandas as pd
import os
import pickle

def find_vertex_with_max_edge_sum(cellmap, vertex_list=None):
    """
    Identifies the vertex with the highest sum of connected edge lengths.
    Optionally restricts the search to vertices in `vertex_list`.
    
    Parameters:
    - cellmap: The cell map containing vertex and edge data.
    - vertex_list: Optional list of vertex IDs to limit the search (default is None, which includes all vertices).
    
    Returns:
    - The vertex ID connected to the highest sum of edge lengths within the specified list, if provided.
    """
    max_sum = -np.inf
    chosen_vertex = None

    # If vertex_list is provided, filter to those vertices only
    vertices_to_check = cellmap.vert_df if vertex_list is None else cellmap.vert_df.loc[vertex_list]

    # Loop through each vertex in the specified set
    for vertex_id, vertex_data in vertices_to_check.iterrows():
        # Get the x and y coordinates of the current vertex
        vertex_x = vertex_data['x']
        vertex_y = vertex_data['y']
        
        # Find all edges that are connected to this vertex
        connected_edges = cellmap.edge_df[
            (cellmap.edge_df['srce'] == vertex_id) | (cellmap.edge_df['trgt'] == vertex_id)
        ]
        
        # Sum the lengths of the connected edges
        edge_sum = connected_edges['length'].sum()

        # Check if this vertex has the highest sum of edge lengths
        if edge_sum > max_sum:
            max_sum = edge_sum
            chosen_vertex = vertex_id

    return chosen_vertex

def collapse_chosen_edge(sheet, chosen_edge):
    """
    Merge the two vertices of `chosen_edge` into a new vertex and
    rewire edges. No plotting, no energy solve, no file I/O.
    Returns: sheet, new_vertex_id
    """
    edge = sheet.edge_df.loc[chosen_edge]
    srce_vertex, trgt_vertex = edge['srce'], edge['trgt']

    # Opposite edge (if any)
    opposites = sheet.edge_df[
        (sheet.edge_df["srce"] == trgt_vertex) & (sheet.edge_df["trgt"] == srce_vertex)
    ]

    # Edges incident to either vertex
    connected_edges_srce = sheet.edge_df[(sheet.edge_df["srce"] == srce_vertex) | (sheet.edge_df["trgt"] == srce_vertex)]
    connected_edges_trgt = sheet.edge_df[(sheet.edge_df["srce"] == trgt_vertex) | (sheet.edge_df["trgt"] == trgt_vertex)]

    # New merged vertex at midpoint
    new_vertex_id = int(max(sheet.vert_df.index)) + 1
    sheet.vert_df.loc[new_vertex_id] = {
        'x': (sheet.vert_df.loc[srce_vertex, 'x'] + sheet.vert_df.loc[trgt_vertex, 'x']) / 2,
        'y': (sheet.vert_df.loc[srce_vertex, 'y'] + sheet.vert_df.loc[trgt_vertex, 'y']) / 2,
        'viscosity': (sheet.vert_df.loc[srce_vertex, 'viscosity'] + sheet.vert_df.loc[trgt_vertex, 'viscosity']) / 2,
        'is_active': 1.0
    }

    # Rewire edges to the new vertex
    connected_edges = pd.concat([connected_edges_srce, connected_edges_trgt]).drop_duplicates()
    for _, e in connected_edges.iterrows():
        if e['srce'] in (srce_vertex, trgt_vertex):
            sheet.edge_df.at[e.name, 'srce'] = new_vertex_id
        if e['trgt'] in (srce_vertex, trgt_vertex):
            sheet.edge_df.at[e.name, 'trgt'] = new_vertex_id

    # Drop the chosen edge (and the opposite if present)
    if not opposites.empty:
        sheet.edge_df.drop([chosen_edge, opposites.index[0]], inplace=True)
    else:
        sheet.edge_df.drop(chosen_edge, inplace=True)

    # Remove old vertices and tidy topology
    sheet.vert_df.drop([srce_vertex, trgt_vertex], inplace=True)
    sheet.reset_index()
    sheet.reset_topo()

    return sheet, new_vertex_id

def split_vertex(cellmap, chosen_vertex, geom, energyContributions_model, distance, retry_attempts=3):
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
            # --- Work on temporary copies ---
            temp_vert_df = cellmap.vert_df.copy()
            temp_edge_df = cellmap.edge_df.copy()

            # --- 1) Find connected edges
            connected_edges = temp_edge_df[
                (temp_edge_df['srce'] == chosen_vertex) |
                (temp_edge_df['trgt'] == chosen_vertex)
            ].copy()

            if connected_edges.empty:
                raise ValueError(f"Vertex {chosen_vertex} has no connected edges.")

            # --- 2) Create new vertex nearby
            new_vert_data = temp_vert_df.loc[chosen_vertex].copy()
            angle = np.random.uniform(0, 2*np.pi)
            dx = distance * np.cos(angle)
            dy = distance * np.sin(angle)
            new_vert_data[cellmap.coords] = temp_vert_df.loc[chosen_vertex, cellmap.coords] + [dx, dy]

            new_vert_index = int(temp_vert_df.index.max()) + 1
            temp_vert_df.loc[new_vert_index] = new_vert_data
    
            source_edge = connected_edges.iloc[0]
        
            # --- 3) Create a pair of edges
            if temp_edge_df.empty:
                template = pd.Series({c: np.nan for c in cellmap.edge_df.columns})
            else:
                template = source_edge.copy().reindex(temp_edge_df.columns, fill_value=np.nan)

            new_edge_index = int(temp_edge_df.index.max()) + 1 if not temp_edge_df.empty else 0
            opposite_edge_index = new_edge_index + 1

            new_edge = template.copy()
            new_edge["srce"], new_edge["trgt"] = chosen_vertex, new_vert_index
            new_edge["face"] = np.nan  # clear before insertion

            opposite_edge = template.copy()
            opposite_edge["srce"], opposite_edge["trgt"] = new_vert_index, chosen_vertex
            opposite_edge["face"] = np.nan  # clear before insertion

            temp_edge_df.loc[new_edge_index] = new_edge
            temp_edge_df.loc[opposite_edge_index] = opposite_edge



            # --- 4) Reassign original edges to closer vertex
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

            # --- 5) Identify "open" faces (whose edge chains don't close)
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

            # filter: only faces touching the new vertex
            open_faces_touching_new = []
            for f in open_faces:
                verts_f = temp_edge_df[temp_edge_df['face'] == f][['srce', 'trgt']].values.ravel()
                if new_vert_index in verts_f or chosen_vertex in verts_f:
                    open_faces_touching_new.append(f)

            print("🔎 Open faces touching new vertex:", open_faces_touching_new)

            if len(open_faces_touching_new) != 2:
                raise ValueError(
                    f"Expected 2 open faces, found {len(open_faces_touching_new)}: {open_faces_touching_new}"
                )

            # --- Find which new edge closes which open face ---
            for f in open_faces_touching_new:
                f_edges = temp_edge_df[temp_edge_df['face'] == f]

                # collect src/trgt vertices
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
                        print(f"⚠️ Face {f}: expected {needed_edge}, "
                              f"but new={new_pair}, opp={opposite_pair}")
                        
            # --- 6) Verify both assigned faces are closed directed cycles; if not, swap once and recheck

            def _face_closes(df, face_id):
                sub = df[df['face'] == face_id][['srce', 'trgt']]
                if sub.empty:
                    return False
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

            # ✅ --- 7) Commit temp results ---
            cellmap.vert_df = temp_vert_df
            cellmap.edge_df = temp_edge_df

            # Now safe to update geometry
            geom.update_all(cellmap)
            cellmap.reset_topo()
            cellmap.reset_index()

            print(f"✅ Successfully divided vertex {chosen_vertex} → new vertex {new_vert_index}")
            return cellmap, chosen_vertex, new_vert_index, new_edge_index, opposite_edge_index

        except Exception as e:
            print(f"⚠️ split_vertex attempt {attempt+1}/{retry_attempts} failed: {e}")
            # rollback original state
            cellmap.vert_df = original_vert_df.copy()
            cellmap.edge_df = original_edge_df.copy()
            geom.update_all(cellmap)
            cellmap.reset_topo()
            cellmap.reset_index()
            attempt += 1

    print("❌ Failed to divide after multiple attempts.")
    return cellmap, chosen_vertex, None, None, None

def final_state_from_solver(solver):
    """Return the last cellmap from solver.history."""
    last_cm = None
    for _, cm in solver.history:
        last_cm = cm
    return last_cm

def find_vertex_by_edge_sum(cellmap, vertex_list=None, mode="min"):
    """
    Identifies the vertex with the lowest or highest sum of connected edge lengths.
    Optionally restricts the search to vertices in `vertex_list`.

    Parameters:
    - cellmap: The cell map containing vertex and edge data.
    - vertex_list: Optional list of vertex IDs to limit the search (default is None, which includes all vertices).
    - mode: "min" to find the vertex with the smallest sum, "max" for the largest sum.

    Returns:
    - The vertex ID connected to the lowest or highest sum of edge lengths within the specified list.
    """
    if mode not in {"min", "max"}:
        raise ValueError("mode must be 'min' or 'max'")

    chosen_vertex = None
    if mode == "min":
        best_sum = np.inf
    else:
        best_sum = -np.inf

    vertices_to_check = cellmap.vert_df if vertex_list is None else cellmap.vert_df.loc[vertex_list]

    for vertex_id, vertex_data in vertices_to_check.iterrows():
        connected_edges = cellmap.edge_df[
            (cellmap.edge_df['srce'] == vertex_id) | (cellmap.edge_df['trgt'] == vertex_id)
        ]
        edge_sum = connected_edges['length'].sum()

        if (mode == "min" and edge_sum < best_sum) or (mode == "max" and edge_sum > best_sum):
            best_sum = edge_sum
            chosen_vertex = vertex_id

    return chosen_vertex

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

def load_simulation_state(folder_path, filename='final_cellmap_state.pkl'):
    file_path = os.path.join(folder_path, filename)
    
    # Load the cellmap object using pickle
    with open(file_path, 'rb') as f:
        cellmap = pickle.load(f)
    
    print(f"Simulation state loaded from {file_path}")
    return cellmap
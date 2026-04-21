import numpy as np

# Initialize physical parameters and boundary conditions
# Also provides analytical reference pressure and flux fields

def initPhysicalParams(cell_struct, face_struct, Lx, Ly, Lz, bc_option):

    # Physical constants 
    K_tensor = np.eye(3)
    phi_vals = 0.3
    rho_vals = 1000.0
    g_val = 0.0
    gravity_dir = np.array([0.0, 0.0, -1.0])
    tol = 1e-6

    n_faces = len(face_struct)

    # Face centers
    f_centers = np.zeros((n_faces, 3))
    for f in range(n_faces):
        f_centers[f, :] = np.array(face_struct[f]["center"]).reshape(3,)

    # Boundary detection 
    west_idx = np.where(np.abs(f_centers[:, 0] - 0.0) < tol)[0]
    east_idx = np.where(np.abs(f_centers[:, 0] - Lx) < tol)[0]

    south_idx = np.where(np.abs(f_centers[:, 1] - 0.0) < tol)[0]
    north_idx = np.where(np.abs(f_centers[:, 1] - Ly) < tol)[0]

    bottom_idx = np.where(np.abs(f_centers[:, 2] - 0.0) < tol)[0]
    top_idx = np.where(np.abs(f_centers[:, 2] - Lz) < tol)[0]

    BC_Dirichlet_map = {}
    BC_Neumann_map = {}

    if bc_option == "linear":

        grad_pref = np.array([-1.0 / Lx, 0.0, 0.0])
        m_ref_vec = -K_tensor @ grad_pref

        # Dirichlet
        for f in west_idx:
            BC_Dirichlet_map[int(f)] = 1.0
        for f in east_idx:
            BC_Dirichlet_map[int(f)] = 0.0

        # Neumann
        neumann_faces = np.concatenate(
            [south_idx, north_idx, bottom_idx, top_idx]
        )
        for f in neumann_faces:
            BC_Neumann_map[int(f)] = 0.0

    # TO DO: NEED TO DEBUG
    elif bc_option == "corner2corner":

        grad_pref = np.array([-1.0 / Lx, -1.0 / Ly, -1.0 / Lz])
        m_ref_vec = -K_tensor @ grad_pref

        boundary_faces = np.unique(
            np.concatenate(
                [west_idx, east_idx, south_idx, north_idx, bottom_idx, top_idx]
            )
        )

        n_cells = len(cell_struct)
        cell_centers = np.zeros((n_cells, 3))

        for c in range(n_cells):
            cell_centers[c, :] = np.array(cell_struct[c]["center"]).reshape(3,)

        inlet_target = np.array([0.0, 0.0, Lz])
        outlet_target = np.array([Lx, Ly, 0.0])

        inlet_cell = np.argmin(np.linalg.norm(cell_centers - inlet_target, axis=1))
        outlet_cell = np.argmin(np.linalg.norm(cell_centers - outlet_target, axis=1))

        inlet_faces_all = np.array(cell_struct[inlet_cell]["faces"])
        outlet_faces_all = np.array(cell_struct[outlet_cell]["faces"])

        inlet_faces = inlet_faces_all[np.isin(inlet_faces_all, boundary_faces)]
        outlet_faces = outlet_faces_all[np.isin(outlet_faces_all, boundary_faces)]

        dirichlet_faces = np.unique(np.concatenate([inlet_faces, outlet_faces]))

        for f in dirichlet_faces:
            xf, yf, zf = face_struct[int(f)]["center"]
            val = grad_pref[0]*xf + grad_pref[1]*yf + grad_pref[2]*zf + 1.0
            BC_Dirichlet_map[int(f)] = val

        neumann_faces = np.setdiff1d(boundary_faces, dirichlet_faces)
        for f in neumann_faces:
            BC_Neumann_map[int(f)] = 0.0

    # Assign cell properties
    for c in range(len(cell_struct)):
        cell_struct[c]["K"] = K_tensor
        cell_struct[c]["phi"] = phi_vals
        cell_struct[c]["rho"] = rho_vals

    # Assign face properties
    for f in range(n_faces):
        face_struct[f]["gravity"] = g_val * gravity_dir
        face_struct[f]["rho"] = rho_vals
        face_struct[f]["BC_flux"] = None

        if f in BC_Dirichlet_map:
            face_struct[f]["BC_pressure"] = BC_Dirichlet_map[f]

        if f in BC_Neumann_map:
            face_struct[f]["BC_flux"] = BC_Neumann_map[f]

    # Analytical flux
    m_ref_faces = np.zeros(n_faces)

    for f in range(n_faces):
        n_f = np.array(face_struct[f]["normal"])
        n_f = n_f / np.linalg.norm(n_f)
        Af = face_struct[f]["area"]

        m_ref_faces[f] = Af * np.dot(m_ref_vec, n_f)

    phys = {
        "K_tensor": K_tensor,
        "grad_pref": grad_pref,
        "m_ref_vec": m_ref_vec,
        "m_ref_faces": m_ref_faces,
    }

    return cell_struct, face_struct, phys


def projectAnalyticalField(cell_struct, face_struct, phys, a, b, c, d):

    nCells = len(cell_struct)
    nFaces = len(face_struct)

    K_tensor = phys["K_tensor"]
    gradp = np.array([a, b, c])

    p_proj = np.zeros(nCells)

    for k in range(nCells):
        xc = np.array(cell_struct[k]["center"])
        p_proj[k] = a*xc[0] + b*xc[1] + c*xc[2] + d

    m_proj = np.zeros(nFaces)

    for f in range(nFaces):
        n_f = np.array(face_struct[f]["normal"])
        n_f = n_f / np.linalg.norm(n_f)
        A_f = face_struct[f]["area"]

        m_proj[f] = -A_f * np.dot(K_tensor @ gradp, n_f)

    return m_proj, p_proj
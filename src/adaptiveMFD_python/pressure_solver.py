import numpy as np
from scipy.sparse import coo_matrix, bmat, diags
from scipy.sparse.linalg import splu, gmres, LinearOperator
from scikits.umfpack import spsolve
from operators import orth

# Pressure solver for adaptive MFD:
# builds global system (M, B, T), applies BCs, and solves via direct or GMRES

def solve_pressure(cell_struct, face_struct, cellMarking,
                   dt_pressure=1.0, g_c=0.0,
                   eps_solver=1e-11, gmres_niter=50,
                   solver_type="direct"):

    n_cells = len(cell_struct)
    n_faces = len(face_struct)
    dim = len(face_struct[0]["center"])

    face_counts = np.array([len(c["faces"]) for c in cell_struct], dtype=int)
    total_nnz = int(np.sum(face_counts * face_counts))

    rows = np.zeros(total_nnz, dtype=int)
    cols = np.zeros(total_nnz, dtype=int)
    vals = np.zeros(total_nnz, dtype=float)
    idx = 0

    gi_cache = {}
    gj_cache = {}
    for nf in np.unique(face_counts):
        ii, jj = np.meshgrid(np.arange(nf), np.arange(nf), indexing="ij")
        gi_cache[nf] = ii.reshape(-1, order="F")
        gj_cache[nf] = jj.reshape(-1, order="F")

    for cc in range(n_cells):

        face_ids = np.asarray(cell_struct[cc]["faces"], dtype=int)
        cell_nf  = face_ids.size

        Cc = np.asarray(cell_struct[cc]["center"]).reshape(-1)
        K = np.asarray(cell_struct[cc]["K"])
        v = cell_struct[cc]["volume"]
        signs = np.asarray(cell_struct[cc]["faces_orientation"]).reshape(-1)

        Cf_mat = np.array([face_struct[f]["center"] for f in face_ids])
        Nf_mat = np.array([face_struct[f]["normal"] for f in face_ids])
        Af_vec = np.array([face_struct[f]["area"] for f in face_ids])

        C = Cf_mat - Cc
        df_norms = np.sqrt(np.sum(C * C, axis=1))
        signf_vec = np.sign(np.sum((C / df_norms[:, None]) * Nf_mat, axis=1))
        N = Af_vec[:, None] * signf_vec[:, None] * Nf_mat

        if cellMarking[cc] == 0:
            td = np.sum(C * (N @ K), axis=1) / np.sum(C * C, axis=1)
            invT = np.diag(1.0 / np.abs(td))
        else:
            # Simple
            t_loc = 6.0 * np.sum(np.diag(K)) / dim
            Q = orth(N / Af_vec[:, None])
            U = np.eye(cell_nf) - Q @ Q.T
            di = np.diag(1.0 / Af_vec)
            invT_reg = (v / t_loc) * (di @ U @ di)
            invT = (C @ np.linalg.solve(K, C.T)) / v + invT_reg

        sign_mat = np.outer(signs, signs)

        gi = face_ids[gi_cache[cell_nf]]
        gj = face_ids[gj_cache[cell_nf]]
        n2 = cell_nf * cell_nf

        rows[idx:idx+n2] = gi.reshape(-1, order="F")
        cols[idx:idx+n2] = gj.reshape(-1, order="F")
        vals[idx:idx+n2] = (sign_mat * invT).reshape(-1, order="F")
        idx += n2

    M = coo_matrix((vals, (rows, cols)), shape=(n_faces, n_faces)).tocsr()
    B = buildBmatrix(cell_struct, face_struct)
    T = buildTmatrix(cell_struct)

    A_full = bmat([
        [M,   -B.T],
        [B,   diags(np.zeros(n_cells))]
    ], format="csr")

    matrix = A_full.copy()

    rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct)
    matrix, rhs_Neumann, _ = neumannBoundary(matrix, face_struct)

    p_n = np.zeros(n_cells)

    g_vec = np.array([0.0, 0.0, -g_c])
    face_normals = np.array([f["normal"] for f in face_struct])
    face_areas = np.array([f["area"] for f in face_struct])
    face_rho = np.array([f["rho"] for f in face_struct])

    f_g = -face_rho * (face_normals @ g_vec) * face_areas

    BC_face_flux_ids = np.array(
        [i for i, s in enumerate(face_struct) if s.get("BC_flux") is not None],
        dtype=int
    )
    BC_face_flux_vals = np.array([face_struct[i]["BC_flux"] for i in BC_face_flux_ids])

    RHS = np.concatenate([f_g + rhs_Dirichlet, (1.0 / dt_pressure) * (T @ p_n)])

    matrix, RHS = enforcePrescribedDOFsStrong(
        BC_face_flux_ids, BC_face_flux_vals, matrix, RHS
    )

    # ===== Solver =====
    # TO DO: Speed up GMRES solver 
    # option 1: direct solver (direct LU solve + iterative refinement)
    if solver_type == "direct":
        sol3 = spsolve(matrix, -RHS)

        for _ in range(3):
            r = matrix @ sol3 + RHS
            sol3 -= spsolve(matrix, r)
    
    # option 2: GMRES
    else:
        num_m_dofs = n_faces
        num_p_dofs = n_cells

        m_dofs = np.arange(num_m_dofs)
        p_dofs = np.arange(num_m_dofs, num_m_dofs + num_p_dofs)

        A_mm = matrix[m_dofs[:, None], m_dofs]
        A_pm = matrix[p_dofs[:, None], m_dofs]

        F_mm = splu(A_mm.tocsc())

        def prec_apply(v):
            r1 = v[:num_m_dofs]
            r2 = v[num_m_dofs:]
            y1 = F_mm.solve(r1)
            y2 = r2 - A_pm @ y1
            return np.concatenate([y1, y2])

        M_prec = LinearOperator(matrix.shape, matvec=prec_apply)

        sol3, info = gmres(
            matrix,
            -RHS,
            M=M_prec,
            rtol=eps_solver,
            atol=0.0,
            restart=30,
            maxiter=gmres_niter
        )

        if info != 0:
            print(f"GMRES info = {info}")

        for _ in range(2):
            r = matrix @ sol3 + RHS
            corr, info_ref = gmres(
                matrix,
                -r,
                M=M_prec,
                x0=np.zeros_like(sol3),
                rtol=eps_solver,
                atol=0.0,
                restart=30,
                maxiter=gmres_niter
            )
            sol3 += corr

    m_num = sol3[:n_faces]
    p_num = sol3[n_faces:]

    return m_num, p_num


def buildBmatrix(cell_struct, face_struct):
    n_cells = len(cell_struct)
    n_faces = len(face_struct)

    total_nnz = sum(len(cell["faces"]) for cell in cell_struct)

    rows = np.zeros(total_nnz, dtype=int)
    cols = np.zeros(total_nnz, dtype=int)
    vals = np.zeros(total_nnz, dtype=float)

    idx = 0
    for k, cell in enumerate(cell_struct):
        face_ids = np.array(cell["faces"]).astype(int)
        signs = np.array(cell["faces_orientation"]).astype(float)

        m = len(face_ids)
        rows[idx:idx+m] = k
        cols[idx:idx+m] = face_ids
        vals[idx:idx+m] = signs
        idx += m

    return coo_matrix((vals, (rows, cols)), shape=(n_cells, n_faces)).tocsr()


def buildTmatrix(cell_struct):
    n_cells = len(cell_struct)
    t_vec = np.zeros(n_cells)

    for i in range(n_cells):
        phi = cell_struct[i]["phi"]
        V = cell_struct[i]["volume"]
        t_vec[i] = phi * V

    return diags(t_vec, 0, shape=(n_cells, n_cells))


def dirichletBoundary(cell_struct, face_struct):
    n_cells = len(cell_struct)
    n_faces = len(face_struct)

    rhs_Dirichlet = np.zeros(n_faces)

    for k in range(n_cells):
        face_ids = np.array(cell_struct[k]["faces"]).astype(int)
        signs = np.array(cell_struct[k]["faces_orientation"]).astype(float)

        for j in range(len(face_ids)):
            f = face_ids[j]
            if face_struct[f].get("BC_pressure") is None:
                continue

            sigma = signs[j]
            p_D = face_struct[f]["BC_pressure"]
            rhs_Dirichlet[f] = sigma * p_D

    return rhs_Dirichlet


def neumannBoundary(A, face_struct):
    n_faces = len(face_struct)
    rhs_BC = np.zeros(n_faces)

    f_ids = np.array([i for i, x in enumerate(face_struct) if x.get("BC_flux") is not None], dtype=int)
    f_vals = np.array([face_struct[i]["BC_flux"] for i in f_ids])

    A = A.tolil()
    A[f_ids, :] = 0
    A[:, f_ids] = 0
    A[f_ids, f_ids] = 1.0
    A = A.tocsr()

    rhs_BC[f_ids] = f_vals
    return A, rhs_BC, f_ids


def buildGravityRHS(face_struct, g):
    n_faces = len(face_struct)
    f_g = np.zeros(n_faces)

    for f in range(n_faces):
        rho = face_struct[f]["rho"]
        n = np.array(face_struct[f]["normal"]).reshape(-1)
        area_f = face_struct[f]["area"]
        f_g[f] = -rho * np.dot(np.array([0.0, 0.0, -g]), n) * area_f

    return f_g


def enforcePrescribedDOFsStrong(prescribedIdx, prescribedVal, A, b):
    nUnknowns = A.shape[0]

    if np.isscalar(prescribedVal):
        prescribedVal = np.full(len(prescribedIdx), prescribedVal)
    else:
        prescribedVal = np.array(prescribedVal).reshape(-1)

    isFree = np.ones(nUnknowns, dtype=bool)
    isFree[prescribedIdx] = False

    selectFree = diags(isFree.astype(float), 0, shape=(nUnknowns, nUnknowns))
    selectPrescribed = diags((~isFree).astype(float), 0, shape=(nUnknowns, nUnknowns))

    xPrescribed = np.zeros(nUnknowns)
    xPrescribed[prescribedIdx] = prescribedVal

    A_freeRows = selectFree @ A
    b = b - A_freeRows @ xPrescribed

    b[prescribedIdx] = prescribedVal

    A = A_freeRows @ selectFree + selectPrescribed
    return A.tocsr(), b


def block_prec(r, F_mm, A_pm, F_S, num_m_dofs):
    r1 = r[:num_m_dofs]
    r2 = r[num_m_dofs:]

    y1 = F_mm.solve(r1)
    y2 = F_S.solve(r2 - A_pm @ y1)

    return np.concatenate([y1, y2])


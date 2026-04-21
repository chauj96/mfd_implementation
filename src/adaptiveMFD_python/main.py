import numpy as np
from mesh_loader import load_mesh
from physics import initPhysicalParams, projectAnalyticalField
from operators import createMmatrix, createBmatrix
from classification import classify_cells
from pressure_solver import solve_pressure
from io_utils import print_pressure_err, plot_pressure_err

# ===== Step 1: Load mesh =====
cell_struct, face_struct, vertices, Lx, Ly, Lz = load_mesh("twoFaults")

# Set analytical linear pressure field ( f(x,y,z) = ax + by + cz + d )
a = -1.0 / Lx
b = 0.0
c = 0.0
d = 1.0

# ===== Step 2: Physical/discrete operator setup =====
g_c = 0.0
dt_pressure = 1.0
 
cell_struct, face_struct, phys = initPhysicalParams(cell_struct, face_struct, Lx, Ly, Lz, bc_option="linear")
cell_struct = createMmatrix(cell_struct, face_struct, ip_type="tpfa")
cell_struct = createBmatrix(cell_struct)
m_proj, p_proj = projectAnalyticalField(cell_struct, face_struct, phys, a, b, c, d)


# ===== Step 3: Classify cells and solve a pressure field =====
# Solver setup
tol_list = np.array([1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
n_tol = len(tol_list)
eps_solver = 1e-11
gmres_niter = 200
solver_type = "direct"  # "direct" or "iterative"

results = []
for tol in tol_list:

    cellMarking = classify_cells(cell_struct, face_struct, m_proj, p_proj, vertices, a, b, c, d, tol)
    m_num, p_num = solve_pressure(cell_struct, face_struct, cellMarking, dt_pressure, g_c, eps_solver, gmres_niter, solver_type)

    flux_rel_err = np.linalg.norm(m_num - m_proj) / np.linalg.norm(m_proj)
    flux_abs_err = np.linalg.norm(m_num - m_proj)

    results.append([tol, flux_rel_err, flux_abs_err])

# Check flux relative/absolute error
print_pressure_err(results)
plot_pressure_err(results)
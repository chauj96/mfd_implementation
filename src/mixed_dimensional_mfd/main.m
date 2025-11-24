% Fracture is placed right middle of the domain (vertical case)
% Mixed boundary: Dirichlet (left, right) and Neumann (top, bottom)
% No gravity

% clear; clc; close all;

case_type = 'structured';
ip_type = 'tpfa';
dt = 1;
rho = 1000;
g_c = 0.0;
Lx = 1.0; Lz = 1.0;
nx = 4; nz = 2;  % can change the number of cells
n_cells = nx * nz;


% Build struct for 2D cells
if strcmp(case_type, 'structured')
    [cell_struct_2d, face_struct_2d, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz);
else
    domain = @MbbDomain;
    n_cells = 1e+4;
    [cell_struct_2d, face_struct_2d, vertices, cells] = buildPolyGrid(domain, n_cells);
end


[cell_struct_2d, face_struct_2d] = buildCellStruct(cell_struct_2d, face_struct_2d, Lx, nx);
[cell_struct_2d, face_struct_2d, phys] = initPhysicalParams(cell_struct_2d, face_struct_2d, Lx, Lz, case_type);

% Build struct for 1D fractures
cell_struct_1d = build1DStruct(vertices, cell_struct_2d, nz, Lx, Lz);


%% 2D 
% Assemble matrices of 2D geometry
M_2d = buildMmatrixParametric(cell_struct_2d, face_struct_2d, ip_type);
B_2d = buildBmatrix(cell_struct_2d, face_struct_2d);
A_2d = [M_2d, -B_2d';B_2d, zeros(length(cell_struct_2d), length(cell_struct_2d))];

% Apply boundary condition 
rhs_dirichlet = dirichletBoundary(cell_struct_2d, face_struct_2d);
[A_2d, rhs_neumann] = neumannBoundary(A_2d, cell_struct_2d, face_struct_2d);
rhs_BC_2d = [rhs_dirichlet + rhs_neumann; zeros(length(cell_struct_2d),1)];

%% 1D
M_1d = buildMmatrixParametric_1D(cell_struct_1d, vertices, ip_type);
B_1d = buildBmatrix_1d(cell_struct_1d);
A_1d = [M_1d, -B_1d';B_1d, zeros(length(cell_struct_1d), length(cell_struct_1d))];

% Apply boundary condition (for project case, there is no dirichelt on fracture)
[A_1d, rhs_neumann_1d] = neumannBoundary_1d(A_1d, cell_struct_1d);
rhs_BC_1d = [rhs_neumann_1d; zeros(length(cell_struct_1d),1)];

%% Assemble 2D + 1D
M_total = [A_2d, zeros(size(A_2d,1), size(A_1d,2)) ; zeros(size(A_1d,1), size(A_2d,2)), A_1d];

% Apply coupling terms (add T and T^t matrices, Neumann coupling to M_2d)
fracture_aperture = 0;
M_total = applyCoupling(M_total, cell_struct_2d, cell_struct_1d, fracture_aperture);

% Solution structure: [m_2d, p_2d, m_1d, p_1d]
% We can compare against to analytical solution (mass flux): phys.m_ref_faces 
sol = M_total \ -[rhs_BC_2d;rhs_BC_1d]
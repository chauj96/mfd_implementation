clear; clc;

addpath('PolyMesher/');

% case options: 'structured' or 'unstructured'
case_type = 'unstructured';

dt = 1;
nx = 50;
nz = 50;
Lx = 2.0;
Lz = 2.0;
rho = 1000;
g_c = 1e-6;
% g_c = 0;

if strcmp(case_type, 'structured')
    [cell_struct, face_struct, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz);

elseif strcmp(case_type, 'unstructured')
    domain = @MbbDomain; % this domain also sets with Lx = 2, Lz = 2 -> BdBox = [0 2 0 2];
    n_cells = 1000; % we can change the number of cells

    [cell_struct, face_struct, vertices, cells] = buildPolyGrid(domain, n_cells);
end

% Step 1: Initialize physical properties
[cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct, Lx, Lz, case_type);

% Step 2: Build matrix M, B, T / Assemble the matrices
M = buildMmatrix(cell_struct, face_struct);
B = buildBmatrix(cell_struct, face_struct);
T = buildTmatrix(cell_struct);

A = [M, -B';B, (1/dt)*T];

% Step 3: Apply boundary conditions (Dirichlet, Neumann)
rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
[A, rhs_Neumann] = neumannBoundary(A, cell_struct, face_struct);

% Step 4: Build RHS
% p_n = ones(length(cell_struct),1) * 1e5; % constant initial pressure
z_top = max(arrayfun(@(c) c.center(2), cell_struct)); 

% Hydrostatic pressure profile
p_n = zeros(length(cell_struct),1);
% for i = 1:length(cell_struct)
%     z_i = cell_struct(i).center(2);
%     p_n(i) = 1e5 + rho * g * (z_top - z_i);
% end

f_g = buildGravityRHS(face_struct, g_c);
rhs = [f_g; (1/dt) * (T * p_n)];

rhs_BC = [rhs_Dirichlet + rhs_Neumann ; zeros(length(cell_struct),1)];


% Step 5: Solve linear system
sol = A \ -(rhs + rhs_BC);

n_faces = length(face_struct);
m_sol = sol(1:n_faces);
p_sol = sol(n_faces+1:end);

% Step 6: Plot the pressure field
plotPressurePolygonal(vertices, cells, p_sol);


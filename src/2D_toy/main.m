clear; clc; close all;

dt = 1;
nx = 1;
nz = 2;
Lx = 1.0;
Lz = 2.0;
rho = 1000;
% g_c = 10.0 * (1.0e-6);
g_c = 0;

[cell_struct, face_struct] = buildStructureGrid(nx, nz, Lx, Lz);

% Step 1: Initialize physical properties
[cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct, Lx, Lz);

% Step 2: Compute transmissibility for each face
face_struct = computeTransmissibility(cell_struct, face_struct);

M_old = buildMmatrix(face_struct);
M = buildMmatrixMFEM(cell_struct, face_struct);

B = buildBmatrix(cell_struct, face_struct);
T = buildTmatrix(cell_struct);
A = [M, -B';B, (1/dt)*T];

rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
[A, rhs_Neumann] = neumannBoundary(A, cell_struct, face_struct);

% Build RHS
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

n_faces = length(face_struct);

% Step 4: Solve linear system
sol = A \ -(rhs + rhs_BC);

m_sol = sol(1:n_faces);
p_sol = sol(n_faces+1:end);
plotPressure(cell_struct, p_sol);


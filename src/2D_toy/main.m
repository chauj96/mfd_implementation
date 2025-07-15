clear; clc; close all;

dt = 1;
nx = 2;
nz = 2;
Lx = 1.0;
Lz = 2.0;
rho = 1000;
% g = 9.8;
g = 0.0;

[cell_struct, face_struct] = buildStructureGrid(nx, nz, Lx, Lz);

% Step 1: Initialize physical properties
[cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct);

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
for i = 1:length(cell_struct)
    z_i = cell_struct(i).center(2);
    p_n(i) = 1e5 + rho * g * (z_top - z_i);
end

f_g = buildGravityRHS(face_struct, 0);
rhs = [f_g; (1/dt) * (T * p_n)];

% Apply boundary condition
% [bc_cells_top, ~] = getCellsAtBoundary(cell_struct, 'top');
% [bc_cells_bottom, ~] = getCellsAtBoundary(cell_struct, 'bottom');

% p_top = 4;
% p_bottom = p_top + rho * g * Lz;
% p_bottom = 2.0;

% bc_cells = [bc_cells_top; bc_cells_bottom];
% bc_values = [repmat(p_top, length(bc_cells_top), 1); repmat(p_bottom, length(bc_cells_bottom), 1)];
% 
% A_mod = A;
% rhs_mod = rhs;
rhs_BC = [rhs_Dirichlet + rhs_Neumann ; zeros(length(cell_struct),1)];

n_faces = length(face_struct);
% for k = 1:length(bc_cells)
%     c = bc_cells(k);
%     val = bc_values(k);
% 
%     row = n_faces + c;
% 
%     A_mod(row, :) = 0;
%     A_mod(:, row) = 0;
%     A_mod(row, row) = 1;
% 
%     rhs_mod(row) = val;
% end


% Step 4: Solve linear system
sol = A \ -(rhs + rhs_BC);

m_sol = sol(1:n_faces);
p_sol = sol(n_faces+1:end);
plotPressure(cell_struct, p_sol);


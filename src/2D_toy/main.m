clear all; clc; close all;

addpath('PolyMesher/');

% case options: 'structured' or 'unstructured'
case_type = 'structured';

dt = 1;
nx = 21;
nz = 21;
Lx = 1.0;
Lz = 1.0;
rho = 1000;
g_c = 0.0 * (1.0e-6);
% g_c = 0;

if strcmp(case_type, 'structured')
    [cell_struct, face_struct, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz);

elseif strcmp(case_type, 'unstructured')
    domain = @MbbDomain; % this domain also sets with Lx = 1, Lz = 1 -> BdBox = [0 1 0 1];
    n_cells = nz*nx; % we can change the number of cells
    [cell_struct, face_struct, vertices, cells] = buildPolyGrid(domain, n_cells);
end

% Step 1: Initialize physical properties
[cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct, Lx, Lz, case_type);

% Step 2: Build matrix M, B, T / Assemble the matrices
% TPFA case
M = buildMmatrix(cell_struct, face_struct, 'tpfa');

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

m_flux = M \ (-B' * p_sol - rhs_Dirichlet);

% Step 6: Plot the pressure field
plotPressurePolygonal(vertices, cells, p_sol, "direct tpfa");

% Step 7: Plot stream lines from flux
%plotStreamlinesFromFlux(cell_struct, face_struct, m_flux);

% (Option) Check perturbation / mesh structure
% figure;
% for c = 1:length(cells)
%     poly = vertices(cells{c}, :);
%     patch(poly(:,1), poly(:,2), 'w', 'EdgeColor', 'k');
% end
% axis equal;
% title('Check perturbation and mesh structure');


clear all; clc; close all;

addpath('PolyMesher/');

% case options: 'structured' or 'unstructured'
case_type = 'structured';

dt = 1;
nx = 81;
nz = 81;
Lx = 1.0;
Lz = 1.0;
rho = 1000;
g_c = 0.0 * (1.0e-6); % or non-zero gravity if desired

if strcmp(case_type, 'structured')
    [cell_struct, face_struct, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz);
elseif strcmp(case_type, 'unstructured')
    domain = @MbbDomain;
    n_cells = nz * nx;
    [cell_struct, face_struct, vertices, cells] = buildPolyGrid(domain, n_cells);
end

% Initialize physical properties
[cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct, Lx, Lz, case_type);

% Build fixed matrices
B = buildBmatrix(cell_struct, face_struct);
T = buildTmatrix(cell_struct);

% Initial pressure
p_n = zeros(length(cell_struct), 1);

% Gravity force RHS
f_g = buildGravityRHS(face_struct, g_c);
rhs_partial = [f_g; (1/dt) * (T * p_n)];
rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);

% Container to store pressure solutions
p_solutions = struct();

% List of inner product types
ip_types = {'tpfa', 'simple', 'general_parametric'};
n_faces = length(face_struct);

solve_full_LS = true;
for i = 1:length(ip_types)
    ip_type = ip_types{i};
    
    % Build M matrix for the selected method
    M = buildMmatrixParametric(cell_struct, face_struct, ip_type);

    % Assemble full system matrix
    A = [M, -B'; B, (1/dt)*T];
    [A, rhs_Neumann] = neumannBoundary(A, cell_struct, face_struct);
    rhs_BC = [rhs_Dirichlet + rhs_Neumann; zeros(length(cell_struct), 1)];
    rhs = rhs_partial + rhs_BC;

    % Solve linear system
    if solve_full_LS

       sol = A \ -rhs;
       m_sol = sol(1:n_faces);
       p_sol = sol(n_faces+1:end);
       
       [m_proj, p_proj] = projectExactSolution(cell_struct, face_struct);

       M_star = A(1:n_faces, 1:n_faces);
       B_star = A(1:n_faces, n_faces+1:end);
       m_flux = -M_star \ (B_star * p_proj + rhs_Dirichlet);
       diff_flux = (m_proj - m_sol)/norm(m_proj);
       diff_flux_norm = norm(diff_flux);
       diff_p = (p_sol - p_proj)/norm(p_proj);
       diff_p_norm = norm(diff_p);
       darcy_residual = M_star * m_proj + (B_star * p_proj + rhs_Dirichlet);
       aka = 0;
    else
        secondary_idx = 1:n_faces;
        primary_idx = n_faces+1:length(rhs_partial);
        sol = schur_solve(A, -rhs, primary_idx, secondary_idx);      
    end
    
    p_sol = sol(n_faces+1:end);
    
    % Store pressure
    p_solutions.(ip_type) = p_sol;
end

% Compute pairwise differences
fprintf('Pairwise norm differences (L2 norm of pressure):\n\n');

for i = 1:length(ip_types)
    for j = i+1:length(ip_types)
        ref = ip_types{i};
        comp = ip_types{j};
        diff_norm = norm(p_solutions.(ref) - p_solutions.(comp)) / norm(p_solutions.(ref));
        fprintf('|| %s - %s || = %.3e\n', ref, comp, diff_norm);
    end
end

for i = 1:length(ip_types)
    plotPressurePolygonal(vertices, cells, p_solutions.(ip_types{i}),ip_types{i});
    title(ip_types{i});
end

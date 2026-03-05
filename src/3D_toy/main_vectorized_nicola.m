


clear all; clc;
addpath(genpath('FACTORIZE'))

% You may select one of the predefined options below.
% Any custom 2D mesh can also be used, provided that it defines:
% vertices, cell-to-vertex connectivity and follows the face-based
% structure required by the solver. The subsequent 3D extrusion and
% discretization steps are independent of how the 2D mesh is generated.

%% Step 1: build 2D meshes

%(Option 1): SPE11B Model (read CSV geometry files)
% [cell2D, face2D, V2, cells2D] = buildMeshFromCSV('csv_files/points_spe11b.csv', 'csv_files/polygons_spe11b.csv');
% H = 100.0; % extrusion height
% Lx = 8400;
% Ly = 1200.02;
% Lz = 100.0;

% (Option 2): Artificial distorted structured mesh for TPFA/MFD tests
H = 1.0;
Lx = 1.0;
Ly = 1.0;
Lz = 1.0;
alpha = 1;
nx = alpha*30;
ny = alpha*30;
[cell2D, face2D, V2, cells2D] = buildTestMesh(nx, ny, Lx, Ly);

%% Step 2: Extrude 2D mesh along vertical direction (H = height) to get 3D mesh
[cell_struct, face_struct, V3, cells3D] = extend2Dto3D(V2, cells2D, H);

%% Step 3: Assign physical values and get projection of analytical solution (flux and pressure)
ip_type = 'tpfa';
tol = 1e-9; % Tolerance for residual-based cell classification
g_c = 0.0;
dt = 1;

% Analytical solution (Dirichlet: pL = 1, pR = 0 / Neumann - no flow: rest of faces)
a = -1/Lx;
b = 0;
c = 0;
d = 1;

[cell_struct, face_struct, phys] = initPhysicalParams3D(cell_struct, face_struct, Lx, Ly, Lz);
[m_proj, p_proj] = projectAnalyticalField3D(cell_struct, face_struct, phys, a, b, c, d); % Get projection of flux and pressure

%% Step 4: Apply cell classification
cell_struct = createMmatrix(cell_struct, face_struct, ip_type);
cell_struct = createBmatrix(cell_struct);

% Precompute d_K for all faces at once (vectorized, avoids per-cell call)
face_centers = reshape([face_struct.center], 3, [])'; % n_faces x 3
d_all = a*face_centers(:,1) + b*face_centers(:,2) + c*face_centers(:,3) + d; % n_faces x 1

res_3D   = zeros(length(face_struct), 1);
res_count = zeros(length(face_struct), 1); % track contributions per face

flux_scale = norm(m_proj);

for cn = 1:length(cell_struct)
    face_ids = cell_struct(cn).faces;
    signs    = cell_struct(cn).faces_orientation(:);

    M_K = cell_struct(cn).M;
    B_K = cell_struct(cn).B;

    mK  = m_proj(face_ids);
    pK  = p_proj(cn);
    d_K = signs .* d_all(face_ids); % vectorized local Dirichlet (no inner loop)

    R_K = M_K * mK - B_K * pK + d_K; % residual equation
    res_3D(face_ids)    = res_3D(face_ids)    + R_K;
    res_count(face_ids) = res_count(face_ids) + 1;
end

% Average contributions at shared faces (avoids inflation with many facets)
res_3D = res_3D ./ flux_scale; max(res_count, 1);

% Vectorized cell marking via accumarray: any face of the cell exceeds tol → mark 1
face_exceeds = abs(res_3D) > tol; % logical flag per face
% Build cell-to-face incidence without nested cell arrays
face_counts  = arrayfun(@(c) length(c.faces), cell_struct(:)); % n_cells x 1
all_face_ids = cell2mat(arrayfun(@(c) c.faces(:), cell_struct(:), 'UniformOutput', false));
all_cell_ids = repelem((1:length(cell_struct))', face_counts);  % repeat cell index for each of its faces
cellMarking_3D = accumarray(all_cell_ids, face_exceeds(all_face_ids), [length(cell_struct) 1], @any);
cellMarking_3D = double(cellMarking_3D);

tpfa_count = length(cell_struct) - sum(cellMarking_3D);
fprintf("TPFA cells = %d out of %d total cells\n", tpfa_count, length(cell_struct));

% Export to VTU file and save the file
outDir = 'output';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

filename = fullfile(outDir, 'mesh.vtu');
writeExtrudedMeshVTP(filename, V3, cell_struct, face_struct, cellMarking_3D, 'cellMarking');

%% Step 5: Solve the global system after the classification
n_cells = length(cell_struct);
n_faces = length(face_struct);
dim = length(face_struct(1).center);

total_nnz = sum(arrayfun(@(c) length(c.faces)^2, cell_struct));
rows = zeros(total_nnz, 1);
cols = zeros(total_nnz, 1);
vals = zeros(total_nnz, 1);
idx = 0;

for c = 1:n_cells

    face_ids = cell_struct(c).faces;
    cell_nf = length(face_ids);
    Cc = cell_struct(c).center(:);
    K = cell_struct(c).K;
    v = cell_struct(c).volume;
    signs = cell_struct(c).faces_orientation;

    % build local geometry (vectorized over faces)
    Cf_mat = reshape([face_struct(face_ids).center], dim, cell_nf)'; % cell_nf x dim
    Nf_mat = reshape([face_struct(face_ids).normal], dim, cell_nf)'; % cell_nf x dim
    Af_vec = [face_struct(face_ids).area]';                          % cell_nf x 1

    C = Cf_mat - Cc';                                  % cell_nf x dim
    df_norms = sqrt(sum(C.^2, 2));                     % cell_nf x 1
    signf_vec = sign(sum((C ./ df_norms) .* Nf_mat, 2)); % cell_nf x 1
    assert(all(signf_vec == signs(:)), 'Orientation mismatch in cell %d', c);

    N = Af_vec .* signf_vec .* Nf_mat;   % cell_nf x dim
    a = Af_vec;                           % cell_nf x 1

    % SELECT SCHEME
    if cellMarking_3D(c) == 0
        % TPFA
        td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
        invT = diag(1 ./ abs(td));

    else
        % SIMPLE MFD
        t  = 6 * sum(diag(K)) / dim;
        Q  = orth(N ./ a);
        U  = eye(cell_nf) - Q * Q';
        di = diag(1 ./ a);
        invT_reg = (v / t) * (di * U * di);
        invT = (C * (K \ C')) / v + invT_reg;

        % General parametric MFD
        % W  = N * K * N';
        % Qn = orth(N);
        % P  = eye(cell_nf) - Qn * Qn';
        % diW = diag(1 ./ diag(W));
        % invT_reg = (v / cell_nf) * (P * diW * P);
        % invT = (C * (K \ C')) / v + invT_reg;
    end


    % Global assembly (vectorized outer product, no i/j double loop)
    sign_mat  = signs(:) * signs(:)';          % cell_nf x cell_nf sign matrix
    [gi, gj]  = ndgrid(face_ids, face_ids);    % global row/col index grids
    n2        = cell_nf^2;
    rows(idx+1:idx+n2) = gi(:);
    cols(idx+1:idx+n2) = gj(:);
    vals(idx+1:idx+n2) = reshape(sign_mat .* invT, [], 1);
    idx = idx + n2;
end

M = sparse(rows, cols, vals, n_faces, n_faces);
B = buildBmatrix(cell_struct, face_struct);
T = buildTmatrix(cell_struct);

A_full = [M, -B';
          B, (1/dt)*T];


matrix = A_full;


rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);

t = tic;
[A_full, rhs_Neumann, f_ids] = neumannBoundary(A_full, [], face_struct);
fprintf('neumannBoundary time: %.6f s\n', toc(t));


p_n = zeros(n_cells, 1);
f_g = buildGravityRHS(face_struct, g_c);

rhs = [f_g; (1/dt)*(T*p_n)];
rhs_total = rhs + [rhs_Dirichlet + rhs_Neumann; zeros(n_cells,1)];





BC_face_flux.ids = find(~arrayfun(@(s) isempty(s.BC_flux), face_struct));
BC_face_flux.vals = [face_struct(BC_face_flux.ids).BC_flux];

RHS = [ f_g + rhs_Dirichlet; ...
        (1/dt)*(T*p_n)];

t = tic;
[matrix, RHS] = enforcePrescribedDOFsStrong( BC_face_flux.ids, ...
                                             BC_face_flux.vals, ...
                                             matrix, ...
                                             RHS);
fprintf('enforcePrescribedDOFsStrong time: %.6f s\n', toc(t));





tic;
num_m_dofs = length(face_struct);
num_p_dofs = length(cell_struct);

m_dofs = 1:num_m_dofs;

p_dofs = num_m_dofs + 1:num_m_dofs + num_p_dofs;

A_mm = matrix(m_dofs, m_dofs);
A_mp = matrix(m_dofs, p_dofs);
A_pm = matrix(p_dofs, m_dofs);
A_pp = matrix(p_dofs, p_dofs);


S_approx = A_pp - A_pm * spdiags(1./diag(A_mm),0,num_m_dofs,num_m_dofs) * A_mp;

F_mm = factorize(A_mm);
F_S = factorize(S_approx);

%%%%%%%%%%%%%%

function y = block_prec(r,F_mm, A_pm, F_S, num_m_dofs)

r1 = r(1:num_m_dofs);
r2 = r(num_m_dofs+1:end);

y1 = F_mm \ r1;
y2 = F_S \ (r2 - A_pm*y1);

y = [y1; y2];

%%%%%%%%%%%%%%%%

end
t_setup = toc;

t = tic;
[sol3, flag, total_iters, error] = gmres_r(@(v) matrix * v, -RHS, 1e-6, 100, 1, @(v) block_prec(v, F_mm, A_pm, F_S, num_m_dofs), 0*RHS);
t_solve = toc(t);
fprintf('solve linear system time: %.6f s\n', t_solve);


true_relative_residual_norm = norm(-RHS - matrix*sol3) / norm(RHS);


m_num = sol3(1:n_faces);
p_num = sol3(n_faces+1:end);

% Maximum error / Relative error
[maxErr, maxIdx] = max(abs(m_num - m_proj));
relErr = norm(m_num - m_proj) / norm(m_proj);
fprintf('Relative error = %e\n', relErr);

fprintf('m DOFs     : %10d\n', num_m_dofs);
fprintf('p DOFs     : %10d\n', num_p_dofs);
fprintf('GMRES iters: %10d\n', total_iters);
fprintf('True rel res: %12.3e\n', true_relative_residual_norm);
fprintf('Rel error   : %12.3e\n', relErr);
fprintf('Setup time  : %12.2f s\n', t_setup);
fprintf('Solve time  : %12.2f s\n', t_solve);
% clear all; clc;

% You may select one of the predefined options below.
% Any custom 2D mesh can also be used, provided that it defines:
% vertices, cell-to-vertex connectivity and follows the face-based
% structure required by the solver. The subsequent 3D extrusion and
% discretization steps are independent of how the 2D mesh is generated.

%% Step 1: build 2D meshes 

%(Option 1): SPE11B Model (read CSV geometry files)
[cell2D, face2D, V2, cells2D] = buildMeshFromCSV('csv_files/points_spe11b.csv', 'csv_files/polygons_spe11b.csv');
H = 100.0; % extrusion height
Lx = 8400;
Ly = 1200.02;
Lz = 100.0;

% (Option 2): Artificial distorted structured mesh for TPFA/MFD tests
% H = 1.0;
% Lx = 1.0;
% Ly = 1.0;
% Lz = 1.0;
% nx = 30;
% ny = 30;
% [cell2D, face2D, V2, cells2D] = buildTestMesh(nx, ny, Lx, Ly);

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

rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
[A_full, rhs_Neumann] = neumannBoundary(A_full, cell_struct, face_struct);

p_n = zeros(n_cells, 1);
f_g = buildGravityRHS(face_struct, g_c);

rhs = [f_g; (1/dt)*(T*p_n)];
rhs_total = rhs + [rhs_Dirichlet + rhs_Neumann; zeros(n_cells,1)];


%% Solve with GMRES preconditioned by a Riesz-map block-diagonal preconditioner
%
%  The saddle-point system is:
%       [ M    -B' ] [m]   [f]
%       [ B  (1/dt)T ] [p] = [g]
%
%  Riesz-map (block-diagonal) preconditioner:
%       P = blkdiag( M ,  S )
%  where S = B*M^{-1}*B' + (1/dt)*T  (pressure Schur complement).
% 

rhs_solve = -rhs_total;
LS_type = 'direct';
if strcmp(LS_type, 'iterative')
    % block: M
    M_sym    = (M + M') / 2;
    % d_mean_M = mean(abs(diag(M_sym)));
    % M_sym    = M_sym + 1e-6 * d_mean_M * speye(n_faces);
    
    fprintf('M  nnz = %d,  n_faces = %d\n', nnz(M_sym), n_faces);
    % Conditioning diagnostics via eigenvalue ratio (lambda_max / lambda_min)
    lam_max_M = eigs(M_sym, 1, 'largestabs',  'Tolerance', 1e-4, 'MaxIterations', 300);
    lam_min_M = eigs(M_sym, 1, 'smallestabs', 'Tolerance', 1e-4, 'MaxIterations', 300);
    fprintf('M  cond estimate (lmax/lmin) = %e  [lmax=%e, lmin=%e]\n', ...
            abs(lam_max_M) / max(abs(lam_min_M), eps), lam_max_M, lam_min_M);
    
    % Escalating ichol: retry with increasing diagonal shift until success
    apply_M_inv = [];
    for alpha = [0, 1e-3, 1e-2, 1e-1, 1.0]
        try
            L_M = ichol(M_sym + alpha * d_mean_M * speye(n_faces), ...
                        struct('type','ict','droptol',1e-3,'diagcomp',0));
            apply_M_inv = @(v) L_M' \ (L_M \ v);
            if alpha > 0
                fprintf('ichol(M) succeeded with alpha = %e\n', alpha);
            end
            break;
        catch
        end
    end
    if isempty(apply_M_inv)
        warning('ichol failed for M block, using Jacobi preconditioner.');
        dM = diag(M_sym);  dM(dM==0) = 1;
        apply_M_inv = @(v) v ./ dM;
    end
    
    % block: Schur complement S = B*M^{-1}*B' + (1/dt)*T 
    row_sum_M = abs(M_sym) * ones(n_faces, 1);
    dM_inv    = 1 ./ max(row_sum_M, 1e-14 * d_mean_M);
    D_Minv    = spdiags(dM_inv, 0, n_faces, n_faces);
    S_approx  = B * D_Minv * B' + (1/dt) * T;
    S_approx  = (S_approx + S_approx') / 2;
    d_mean_S  = mean(abs(diag(S_approx)));
    S_approx  = S_approx + 1e-6 * d_mean_S * speye(n_cells);  % adaptive mesh-scale shift
    
    fprintf('S  nnz = %d,  n_cells = %d\n', nnz(S_approx), n_cells);
    % Conditioning diagnostics via eigenvalue ratio (lambda_max / lambda_min)
    lam_max_S = eigs(S_approx, 1, 'largestabs',  'Tolerance', 1e-4, 'MaxIterations', 300);
    lam_min_S = eigs(S_approx, 1, 'smallestabs', 'Tolerance', 1e-4, 'MaxIterations', 300);
    fprintf('S  cond estimate (lmax/lmin) = %e  [lmax=%e, lmin=%e]\n', ...
            abs(lam_max_S) / max(abs(lam_min_S), eps), lam_max_S, lam_min_S);
    
    
    % Escalating ichol for S
    apply_S_inv = [];
    for alpha = [0, 1e-3, 1e-2, 1e-1, 1.0]
        try
            L_S = ichol(S_approx + alpha * d_mean_S * speye(n_cells), ...
                        struct('type','ict','droptol',1e-3,'diagcomp',0));
            apply_S_inv = @(v) L_S' \ (L_S \ v);
            if alpha > 0
                fprintf('ichol(S) succeeded with alpha = %e\n', alpha);
            end
            break;
        catch
        end
    end
    if isempty(apply_S_inv)
        warning('ichol failed for Schur block, using Jacobi preconditioner.');
        dS = diag(S_approx);  dS(dS==0) = 1;
        apply_S_inv = @(v) v ./ dS;
    end
    
    % block-diagonal Riesz-map preconditioner
    apply_P_inv = @(v) [ apply_M_inv(v(1:n_faces)); ...
                         apply_S_inv(v(n_faces+1:end)) ];
    
    % GMRES
    gmres_tol     = tol * 1e-3;
    gmres_maxit   = 500;
    gmres_restart = 75;
    
    [sol, gmres_flag, gmres_relres, gmres_iter] = ...
        gmres(@(v) A_full*v, rhs_solve, gmres_restart, gmres_tol, gmres_maxit, apply_P_inv);
    
    if gmres_flag == 0
        fprintf('GMRES converged: rel. residual = %e, outer/inner iter = %d/%d\n', ...
                gmres_relres, gmres_iter(1), gmres_iter(2));
    else
        warning('GMRES did not converge (flag=%d), rel. residual = %e', ...
                gmres_flag, gmres_relres);
    end

else
    sol = A_full \ (-rhs_total);
end

m_num = sol(1:n_faces);
p_num = sol(n_faces+1:end);

% Maximum error / Relative error
[maxErr, maxIdx] = max(abs(m_num - m_proj));
relErr = norm(m_num - m_proj) / norm(m_proj);
fprintf('Relative error = %e\n', relErr);

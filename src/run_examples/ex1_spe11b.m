clear all; clc;
addpath(genpath('FACTORIZE'))

%% Step 1: build 2D meshes
[cell2D, face2D, V2, cells2D] = buildMeshFromCSV('csv_files/points_spe11b.csv', 'csv_files/polygons_spe11b.csv');
H = 100.0; % extrusion height
Lx = 8400;
Ly = 1200.02;
Lz = 100.0;

%% Step 2: Extrude 2D mesh along vertical direction (H = height) to get 3D mesh
[cell_struct, face_struct, V3, cells3D] = extend2Dto3D(V2, cells2D, H);

%% Step 3: Assign physical values and get projection of analytical solution (flux and pressure)
ip_type = 'tpfa';
solver_type = 'iterative';
tol_values = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
n_tol = length(tol_values);
abs_l2_errors = zeros(n_tol,1);
rel_l2_errors = zeros(n_tol,1);
tpfa_counts = zeros(n_tol,1);

nnz_adaptive_vals = zeros(n_tol,1);
sparsity_reduction_vals = zeros(n_tol,1);
gmres_iters_vals = zeros(n_tol,1);

% Iterative solver parameters
eps_solver = 1.0e-11;
gmres_niter = 1000;
g_c = 0.0;
dt = 1;

% Analytical solution (Dirichlet: pL = 1, pR = 0 / Neumann - no flow: rest of faces)
a = -1/Lx;
b = 0;
c = 0;
d = 1;

[cell_struct, face_struct, phys] = initPhysicalParams3D(cell_struct, face_struct, Lx, Ly, Lz, 'identity', 'linear');
[m_proj, p_proj] = projectAnalyticalField3D(cell_struct, face_struct, phys, a, b, c, d);

cell_struct = createMmatrix(cell_struct, face_struct, ip_type);
cell_struct = createBmatrix(cell_struct);

% Precompute d_K for all faces at once (vectorized, avoids per-cell call)
face_centers = reshape([face_struct.center], 3, [])';
d_all = a*face_centers(:,1) + b*face_centers(:,2) + c*face_centers(:,3) + d;

n_cells = length(cell_struct);
n_faces = length(face_struct);
dim = length(face_struct(1).center);
total_nnz = sum(arrayfun(@(cc) length(cc.faces)^2, cell_struct));

B = buildBmatrix(cell_struct, face_struct);
T = buildTmatrix(cell_struct);
rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);

p_n = zeros(n_cells, 1);
f_g = buildGravityRHS(face_struct, g_c);

BC_face_flux.ids = find(~arrayfun(@(s) isempty(s.BC_flux), face_struct));
BC_face_flux.vals = [face_struct(BC_face_flux.ids).BC_flux];

RHS_base = [f_g + rhs_Dirichlet; ...
            (1/dt)*(T*p_n)];

%% ===== FULL MFD REFERENCE MATRIX / SOLVE =====
fprintf('\n=== Building and solving FULL MFD system ===\n');

rows = zeros(total_nnz, 1);
cols = zeros(total_nnz, 1);
vals = zeros(total_nnz, 1);
idx = 0;

for cc = 1:n_cells
    face_ids = cell_struct(cc).faces;
    cell_nf = length(face_ids);
    Cc = cell_struct(cc).center(:);
    K = cell_struct(cc).K;
    v = cell_struct(cc).volume;
    signs = cell_struct(cc).faces_orientation(:);

    % build local geometry (vectorized over faces)
    Cf_mat = reshape([face_struct(face_ids).center], dim, cell_nf)';
    Nf_mat = reshape([face_struct(face_ids).normal], dim, cell_nf)';
    Af_vec = [face_struct(face_ids).area]';

    C = Cf_mat - Cc';
    df_norms = sqrt(sum(C.^2, 2));
    signf_vec = sign(sum((C ./ df_norms) .* Nf_mat, 2));
    assert(all(signf_vec == signs(:)), 'Orientation mismatch in cell %d', cc);

    N = Af_vec .* signf_vec .* Nf_mat;
    af = Af_vec;

    % ===== FULL MFD INNER PRODUCT =====
    % SIMPLE MFD
    t_loc = 6 * sum(diag(K)) / dim;
    Q  = orth(N ./ af);
    U  = eye(cell_nf) - Q * Q';
    di = diag(1 ./ af);
    invT_reg = (v / t_loc) * (di * U * di);
    invT = (C * (K \ C')) / v + invT_reg;

    % General Parametric
    % W  = N * K * N';
    % Qn = orth(N);
    % P  = eye(cell_nf) - Qn * Qn';
    % diW = diag(1 ./ diag(W));
    % invT_reg = (v / cell_nf) * (P * diW * P);
    % invT = (C * (K \ C')) / v + invT_reg;

    sign_mat = signs * signs';
    [gi, gj] = ndgrid(face_ids, face_ids);
    n2 = cell_nf^2;

    rows(idx+1:idx+n2) = gi(:);
    cols(idx+1:idx+n2) = gj(:);
    vals(idx+1:idx+n2) = reshape(sign_mat .* invT, [], 1);
    idx = idx + n2;
end

M_full = sparse(rows, cols, vals, n_faces, n_faces);
nnz_full = nnz(M_full);
fprintf('nnz(full MFD) = %d\n', nnz_full);

A_full_ref = [M_full, -B';
              B, (1/dt)*T*0];

matrix_full = A_full_ref;
[matrix_full, RHS_full] = enforcePrescribedDOFsStrong(BC_face_flux.ids, ...
                                                       BC_face_flux.vals, ...
                                                       matrix_full, ...
                                                       RHS_base);

num_m_dofs = n_faces;
num_p_dofs = n_cells;
m_dofs = 1:num_m_dofs;
p_dofs = num_m_dofs + 1:num_m_dofs + num_p_dofs;

A_mm = matrix_full(m_dofs, m_dofs);
A_mp = matrix_full(m_dofs, p_dofs);
A_pm = matrix_full(p_dofs, m_dofs);
A_pp = matrix_full(p_dofs, p_dofs);

S_approx = A_pp - A_pm * spdiags(1./diag(A_mm), 0, num_m_dofs, num_m_dofs) * A_mp;

F_mm = factorize(A_mm);
F_S = factorize(S_approx);

t_full = tic;
if strcmp(solver_type, 'direct')
    sol_full = matrix_full \ (-RHS_full);
    total_iters_full = 0;
else
    [sol_full, flag_full, total_iters_full, error_full] = gmres_r( ...
        @(v) matrix_full * v, ...
        -RHS_full, ...
        eps_solver, ...
        gmres_niter, ...
        1, ...
        @(v) block_prec(v, F_mm, A_pm, F_S, num_m_dofs), ...
        0*RHS_full);
end
solve_time_full = toc(t_full);

m_full = sol_full(1:n_faces);
p_full = sol_full(n_faces+1:end);

fprintf('Full MFD GMRES iters: %d\n', total_iters_full);
fprintf('Full MFD solve time  : %.6f s\n', solve_time_full);
fprintf('Full MFD rel error vs projection solution = %.6e\n', norm(m_full - m_proj) / norm(m_proj));

%% ===== ADAPTIVE LOOP =====
for it = 1:n_tol
    tol = tol_values(it);
    fprintf('\n=== Running with tol = %e (iter %d/%d) ===\n', tol, it, n_tol);

    %% Step 4: Apply cell classification
    res_3D = zeros(n_faces, 1);

    for cn = 1:n_cells
        face_ids = cell_struct(cn).faces;
        signs = cell_struct(cn).faces_orientation(:);

        M_K = signs .* cell_struct(cn).M;
        B_K = cell_struct(cn).B;

        mK = signs .* m_proj(face_ids);
        pK = p_proj(cn);
        d_K = signs .* d_all(face_ids);

        DeltaP_K = -B_K * pK + d_K;
        R_K = (M_K * mK - B_K * pK + d_K) / norm(DeltaP_K);
        res_3D(face_ids) = res_3D(face_ids) + R_K;
    end

    % Vectorized cell marking via accumarray
    face_exceeds = abs(res_3D) > tol;
    face_counts = arrayfun(@(cc) length(cc.faces), cell_struct(:));
    all_face_ids = cell2mat(arrayfun(@(cc) cc.faces(:), cell_struct(:), 'UniformOutput', false));
    all_cell_ids = repelem((1:n_cells)', face_counts);

    cellMarking_3D = accumarray(all_cell_ids, face_exceeds(all_face_ids), [n_cells 1], @any);
    cellMarking_3D = double(cellMarking_3D);

    tpfa_count = n_cells - sum(cellMarking_3D);
    fprintf('tol = %.1e | TPFA cells = %d / %d\n', tol, tpfa_count, n_cells);

    % Export classification to VTU
    outDir = 'output_spe11b';
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    filename = fullfile(outDir, sprintf('mesh_l_%d.vtu', it-1));
    writeExtrudedMeshVTP(filename, V3, cell_struct, face_struct, cellMarking_3D, 'cellMarking', 'cell_plot');

    %% Step 5: Solve the global system after classification
    rows = zeros(total_nnz, 1);
    cols = zeros(total_nnz, 1);
    vals = zeros(total_nnz, 1);
    idx = 0;

    for cc = 1:n_cells
        face_ids = cell_struct(cc).faces;
        cell_nf = length(face_ids);
        Cc = cell_struct(cc).center(:);
        K = cell_struct(cc).K;
        v = cell_struct(cc).volume;
        signs = cell_struct(cc).faces_orientation(:);

        % build local geometry (vectorized over faces)
        Cf_mat = reshape([face_struct(face_ids).center], dim, cell_nf)';
        Nf_mat = reshape([face_struct(face_ids).normal], dim, cell_nf)';
        Af_vec = [face_struct(face_ids).area]';

        C = Cf_mat - Cc';
        df_norms = sqrt(sum(C.^2, 2));
        signf_vec = sign(sum((C ./ df_norms) .* Nf_mat, 2));
        assert(all(signf_vec == signs(:)), 'Orientation mismatch in cell %d', cc);

        N = Af_vec .* signf_vec .* Nf_mat;
        af = Af_vec;

        % SELECT SCHEME
        if cellMarking_3D(cc) == 0
            % TPFA
            td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
            invT = diag(1 ./ abs(td));

        else
            % SIMPLE 
            t_loc = 6 * sum(diag(K)) / dim;
            Q  = orth(N ./ af);
            U  = eye(cell_nf) - Q * Q';
            di = diag(1 ./ af);
            invT_reg = (v / t_loc) * (di * U * di);
            invT = (C * (K \ C')) / v + invT_reg;

            % General Parametric
            % W  = N * K * N';
            % Qn = orth(N);
            % P  = eye(cell_nf) - Qn * Qn';
            % diW = diag(1 ./ diag(W));
            % invT_reg = (v / cell_nf) * (P * diW * P);
            % invT = (C * (K \ C')) / v + invT_reg;
        end

        sign_mat = signs * signs';
        [gi, gj] = ndgrid(face_ids, face_ids);
        n2 = cell_nf^2;

        rows(idx+1:idx+n2) = gi(:);
        cols(idx+1:idx+n2) = gj(:);
        vals(idx+1:idx+n2) = reshape(sign_mat .* invT, [], 1);
        idx = idx + n2;
    end

    M = sparse(rows, cols, vals, n_faces, n_faces);

    nnz_adaptive = nnz(M);
    sparsity_reduction = 1 - nnz_adaptive / nnz_full;
    nnz_adaptive_vals(it) = nnz_adaptive;
    sparsity_reduction_vals(it) = sparsity_reduction;

    fprintf('nnz(full MFD) = %d, nnz(adaptive) = %d, reduction = %.2f%%\n', ...
        nnz_full, nnz_adaptive, 100*sparsity_reduction);

    A_full = [M, -B';
              B, (1/dt)*T*0];

    matrix = A_full;

    [matrix, RHS] = enforcePrescribedDOFsStrong(BC_face_flux.ids, ...
                                                BC_face_flux.vals, ...
                                                matrix, ...
                                                RHS_base);

    % Iterative solve for adaptive system
    A_mm = matrix(m_dofs, m_dofs);
    A_mp = matrix(m_dofs, p_dofs);
    A_pm = matrix(p_dofs, m_dofs);
    A_pp = matrix(p_dofs, p_dofs);

    S_approx = A_pp - A_pm * spdiags(1./diag(A_mm), 0, num_m_dofs, num_m_dofs) * A_mp;

    F_mm = factorize(A_mm);
    F_S = factorize(S_approx);

    t_solve = tic;
    if strcmp(solver_type, 'direct')
        sol3 = matrix \ (-RHS);
        total_iters = 0;
    else
        [sol3, flag, total_iters, error] = gmres_r( ...
            @(v) matrix * v, ...
            -RHS, ...
            eps_solver, ...
            gmres_niter, ...
            1, ...
            @(v) block_prec(v, F_mm, A_pm, F_S, num_m_dofs), ...
            0*RHS);
    end
    solve_time = toc(t_solve);
    gmres_iters_vals(it) = total_iters;

    fprintf('Adaptive GMRES iters: %d\n', total_iters);
    fprintf('Adaptive solve time  : %.6f s\n', solve_time);

    m_num = sol3(1:n_faces);
    p_num = sol3(n_faces+1:end);

    % Reference remains analytical projection m_proj, not m_full.
    abs_l2_error = norm(m_num - m_proj);
    rel_l2_error = abs_l2_error / norm(m_proj);

    fprintf('Absolute L2 error = %e\n', abs_l2_error);
    fprintf('Relative L2 error = %e\n', rel_l2_error);

    abs_l2_errors(it) = abs_l2_error;
    rel_l2_errors(it) = rel_l2_error;
    tpfa_counts(it) = tpfa_count;
end

%% ===== FLUX ERROR PLOT =====
figure;

c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];

plot_rel_l2_errors = max(rel_l2_errors, 1e-16);
plot_abs_l2_errors = max(abs_l2_errors, 1e-16);

loglog(tol_values, plot_rel_l2_errors, '-o', ...
    'LineWidth', 2, ...
    'Color', c1, ...
    'MarkerSize', 8, ...
    'DisplayName', 'Relative flux error');
hold on;

loglog(tol_values, plot_abs_l2_errors, '-s', ...
    'LineWidth', 2, ...
    'Color', c2, ...
    'MarkerSize', 8, ...
    'DisplayName', 'Absolute flux error');

loglog(tol_values, tol_values, '--', ...
    'LineWidth', 1.5, ...
    'Color', c1, ...
    'HandleVisibility', 'off');

loglog(tol_values, tol_values, ':', ...
    'LineWidth', 1.5, ...
    'Color', c2, ...
    'HandleVisibility', 'off');

xlabel('Tolerance', 'FontSize', 22);
ylabel('Error', 'FontSize', 22);
title('Mass Flux Error vs Tolerance', 'FontSize', 24);
legend('Location', 'northwest', 'FontSize', 18);

set(gca, ...
    'FontSize', 18, ...
    'LineWidth', 1.5, ...
    'XScale', 'log', ...
    'YScale', 'log');

set(gca, 'Position', [0.13 0.13 0.75 0.75]);
xlim([min(tol_values), max(tol_values)]);
grid on;

set(gcf, ...
    'Units', 'pixels', ...
    'Position', [100 100 1200 900]);

function y = block_prec(r, F_mm, A_pm, F_S, num_m_dofs)

    r1 = r(1:num_m_dofs);
    r2 = r(num_m_dofs+1:end);

    y1 = F_mm \ r1;
    y2 = F_S \ (r2 - A_pm*y1);

    y = [y1; y2];

end
% clear all; clc;
addpath(genpath('FACTORIZE'))

% You may select one of the predefined options below.
% Any custom 2D mesh can also be used, provided that it defines:
% vertices, cell-to-vertex connectivity and follows the face-based
% structure required by the solver. The subsequent 3D extrusion and
% discretization steps are independent of how the 2D mesh is generated.

%% Step 1: build 2D meshes

% (Option 1): SPE11B Model (read CSV geometry files)
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
% alpha = 3;
% nx = alpha*30;
% ny = alpha*30;
% [cell2D, face2D, V2, cells2D] = buildTestMesh(nx, ny, Lx, Ly);

%% Step 2: Extrude 2D mesh along vertical direction (H = height) to get 3D mesh
[cell_struct, face_struct, V3, cells3D] = extend2Dto3D(V2, cells2D, H);

%% Step 3: Assign physical values and get projection of analytical solution (flux and pressure)
ip_type = 'tpfa';
%tol = 1e-4; % Tolerance for residual-based cell classification
% run for a vector of tolerances: 1e-1 .. 1e-9
tol_values = 10.^-(0:9); % [1e-0,1e-1, 1e-2, ..., 1e-9]
n_tol = length(tol_values);
abs_l2_errors = zeros(n_tol,1);
rel_l2_errors = zeros(n_tol,1);
tpfa_counts   = zeros(n_tol,1);
gmres_iters   = zeros(n_tol,1);
cond_M_vals   = zeros(n_tol,1);
nnz_M_vals    = zeros(n_tol,1);
mfd_area_ratio = zeros(n_tol,1);

% We'll loop over tol_values below (start loop after projection)
% per-iteration tolerance (will be set inside the loop)
eps_solver = 1.0e-11;
gmres_niter = 1000;
g_c = 0.0;
dt = 1;

% Analytical solution (Dirichlet: pL = 1, pR = 0 / Neumann - no flow: rest of faces)
a = -1/Lx;
b = 0;
c = 0;
d = 1;

[cell_struct, face_struct, phys] = initPhysicalParams3D(cell_struct, face_struct, Lx, Ly, Lz);
[m_proj, p_proj] = projectAnalyticalField3D(cell_struct, face_struct, phys, a, b, c, d); % Get projection of flux and pressure

 cell_struct = createMmatrix(cell_struct, face_struct, ip_type);
 cell_struct = createBmatrix(cell_struct);

 % Precompute d_K for all faces at once (vectorized, avoids per-cell call)
 face_centers = reshape([face_struct.center], 3, [])'; % n_faces x 3
 d_all = a*face_centers(:,1) + b*face_centers(:,2) + c*face_centers(:,3) + d; % n_faces x 1

% Start loop over tolerances
for it = 1:n_tol
    tol = tol_values(it);
    fprintf('\n=== Running with tol = %e (iter %d/%d) ===\n', tol, it, n_tol);

%% Step 4: Apply cell classification

 res_3D   = zeros(length(face_struct), 1);
 res_count = zeros(length(face_struct), 1); % track contributions per face
 area_vector = [face_struct.area]';
 flux_scale = norm(m_proj);

 for cn = 1:length(cell_struct)
     face_ids = cell_struct(cn).faces;
     signs    = cell_struct(cn).faces_orientation(:);

     M_K = signs .* cell_struct(cn).M;
     B_K = cell_struct(cn).B;

     mK  = signs .* m_proj(face_ids);
     pK  = p_proj(cn);
     d_K = signs .* d_all(face_ids); % vectorized local Dirichlet (no inner loop)

     DeltaP_K = - B_K * pK + d_K;
     % R_K = M_K * mK - B_K * pK + d_K; % residual equation
     R_K = (M_K * mK - B_K * pK + d_K) / norm(DeltaP_K); % scaled residual equation
     res_3D(face_ids)    = res_3D(face_ids)    + R_K;
     % fprintf("Projected local residual norm: %f\n", norm(R_K));
     res_count(face_ids) = res_count(face_ids) + 1;
 end

 % Average contributions at shared faces (avoids inflation with many facets)
 
 fprintf("Projected residual norm: %e\n", norm(res_3D));

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

 % Compute max/min face area ratio for MFD cells
 mfd_cell_ids = find(cellMarking_3D == 1);
 if ~isempty(mfd_cell_ids)
     mfd_face_ids = unique(cell2mat(arrayfun(@(c) cell_struct(c).faces(:), mfd_cell_ids, 'UniformOutput', false)));
     mfd_areas = area_vector(mfd_face_ids);
     mfd_area_ratio(it) = max(mfd_areas) / min(mfd_areas);
 else
     mfd_area_ratio(it) = NaN;
 end
 fprintf('MFD face area ratio (max/min) = %e\n', mfd_area_ratio(it));

 % Export to VTU file and save the file (only for the first tol to avoid many files)
 outDir = 'output';

 if ~exist(outDir, 'dir')
     mkdir(outDir);
 end

 filename = fullfile(outDir, sprintf('mesh_l_%d.vtu', it-1));
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
     af = Af_vec;                          % cell_nf x 1  (renamed to avoid clobbering scalar 'a')

     % SELECT SCHEME
     if cellMarking_3D(c) == 0
         % TPFA
         td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
         invT = diag(1 ./ abs(td));

     else
         % SIMPLE MFD
         t  = 6 * sum(diag(K)) / dim;
         Q  = orth(N ./ af);
         U  = eye(cell_nf) - Q * Q';
         di = diag(1 ./ af);
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
 nnz_M_vals(it) = nnz(M);
 cond_M_vals(it) = condest(M);
 fprintf('nnz(M) = %d,  condest(M) = %e\n', nnz_M_vals(it), cond_M_vals(it));
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

 %%%%%%%%%%%%%%%

 t_setup = toc;

 t = tic;
 [sol3, flag, total_iters, error] = gmres_r(@(v) matrix * v, -RHS, eps_solver, gmres_niter, 1, @(v) block_prec(v, F_mm, A_pm, F_S, num_m_dofs), 0*RHS);
 t_solve = toc(t);
 fprintf('solve linear system time: %.6f s\n', t_solve);

 true_residual_norm = norm(-RHS - matrix*sol3);
 true_relative_residual_norm = norm(-RHS - matrix*sol3) / norm(RHS);


 m_num = sol3(1:n_faces);
 p_num = sol3(n_faces+1:end);

 % Maximum error / Relative error
 [maxErr, maxIdx] = max(abs(m_num - m_proj));



 relErr = norm(m_num - m_proj) / norm(m_proj);
 fprintf('Relative error = %e\n', norm(m_proj));
 fprintf('Error = %e\n', norm(m_num - m_proj));
 fprintf('Relative error = %e\n', relErr);


 % 1. Compute the squared difference
 diff_sq = (m_num - m_proj).^2;

 % 2. Compute the absolute L2 error: sqrt( sum( Area * error^2 ) )
 %abs_l2_error = sqrt(sum(area_vector .* diff_sq));
 abs_l2_error = sqrt(sum(diff_sq));

 % 3. Compute the L2 norm of the reference solution for normalization
 %ref_l2_norm = sqrt(sum(area_vector .* m_proj.^2));
 ref_l2_norm = sqrt(sum(m_proj.^2));

 % 4. Compute Relative Error
 rel_l2_error = abs_l2_error / ref_l2_norm;

 fprintf('Absolute L2 error = %e\n', abs_l2_error);
 fprintf('Relative L2 error = %e\n', rel_l2_error);

    % store errors and solver info for this tolerance
    abs_l2_errors(it) = abs_l2_error;
    rel_l2_errors(it) = rel_l2_error;
    tpfa_counts(it)   = tpfa_count;
    gmres_iters(it)   = total_iters;
    % cond_M_vals and nnz_M_vals already stored above after M assembly

 fprintf('Stored errors for tol=%e (iter %d)\n', tol, it);
 fprintf('GMRES iters: %10d\n', total_iters);
 fprintf('True res: %12.3e\n', true_residual_norm);
 fprintf('True rel res: %12.3e\n', true_relative_residual_norm);
 fprintf('Rel error   : %12.3e\n', relErr);
 fprintf('Setup time  : %12.2f s\n', t_setup);
 fprintf('Solve time  : %12.2f s\n', t_solve);


end % for it

% Plot both errors in one figure with dual y-axes
figure;
left_color  = [0.2 0.4 0.8];   % blue  – relative error axis
right_color = [0.8 0.2 0.2];   % red   – absolute error axis

yyaxis left;
set(gca, 'YColor', left_color);
loglog(tol_values, rel_l2_errors, '-o', 'LineWidth', 1.5, 'Color', left_color,  'DisplayName', 'Relative L2 error');
hold on;
loglog(tol_values, tol_values,    '--',  'LineWidth', 1.2, 'Color', left_color,  'DisplayName', 'y = tol (left)');
set(gca, 'XScale', 'log', 'YScale', 'log');
ylabel('Relative L2 error');

yyaxis right;
set(gca, 'YColor', right_color);
loglog(tol_values, abs_l2_errors, '-s',  'LineWidth', 1.5, 'Color', right_color, 'DisplayName', 'Absolute L2 error');
hold on;
loglog(tol_values, tol_values,    '--',  'LineWidth', 1.2, 'Color', right_color, 'DisplayName', 'y = tol (right)');
set(gca, 'YScale', 'log');
ylabel('Absolute L2 error');

xlabel('tol');
title(sprintf('L2 errors vs tolerance  (||m_{ref}||_{L2} = %.4e)', ref_l2_norm));
legend('Relative L2 error', 'y = tol (left)', 'Absolute L2 error', 'y = tol (right)', 'Location', 'best');
grid on;

% Plot TPFA cell count vs tol
figure;
loglog(tol_values, tpfa_counts, '-^', 'LineWidth', 1.5, 'Color', [0.2 0.6 0.2]);
xlabel('tol');
ylabel('Number of TPFA cells');
title('TPFA cell count vs tolerance');
grid on;
n_total = length(cell_struct);
yline(n_total, '--r', sprintf('Total cells = %d', n_total), 'LabelHorizontalAlignment', 'left');

% Plot: estimated cond(M) vs MFD face area ratio (max/min)
figure;
loglog(mfd_area_ratio, cond_M_vals, '-o', 'LineWidth', 1.5, 'Color', [0.8 0.4 0.0]);
xlabel('Max/min face area ratio (MFD cells)');
ylabel('Estimated cond(M)');
title('Condition number of M vs MFD face area ratio');
grid on;
% annotate each point with the corresponding tol value
for k = 1:n_tol
    if ~isnan(mfd_area_ratio(k))
        text(mfd_area_ratio(k), cond_M_vals(k), sprintf('  tol=%.0e', tol_values(k)), 'FontSize', 7);
    end
end

% Plot: estimated cond(M) [left axis] and GMRES iterations [right axis] vs tol
figure;
yyaxis left;
loglog(tol_values, cond_M_vals, '-o', 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
ylabel('Estimated cond(M)');
set(gca, 'YScale', 'log');

yyaxis right;
semilogx(tol_values, gmres_iters, '-d', 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);
ylabel('GMRES iterations');

xlabel('tol');
title('Condition number of M and GMRES iterations vs tolerance');
legend('condest(M)', 'GMRES iterations', 'Location', 'best');
grid on;

% Plot 2: number of nonzeros in M vs tol
figure;
loglog(tol_values, nnz_M_vals, '-s', 'LineWidth', 1.5, 'Color', [0.5 0.2 0.7]);
xlabel('tol');
ylabel('nnz(M)');
title('Number of nonzeros in M vs tolerance');
grid on;

fprintf('m DOFs     : %10d\n', num_m_dofs);
fprintf('p DOFs     : %10d\n', num_p_dofs);
% Note: preconditioner implemented in separate file `block_prec.m` in the same folder
% to remain compatible with MATLAB versions that don't allow local functions in scripts.

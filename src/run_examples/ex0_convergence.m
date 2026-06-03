clear all; close all; clc;
addpath(genpath('FACTORIZE'))

% ============================================================
%  Select classification approach:
%    true  -> Global  Adaptation (GA): face-accumulated residual
%    false -> Local   Adaptation (LA): cell inf-norm residual
% ============================================================
use_GA = false;

if use_GA
    method_tag   = 'global_adaptation';
    method_label = 'GA';
else
    method_tag   = 'local_adaptation';
    method_label = 'LA';
end

fprintf('============================================================\n');
fprintf('  ADAPTIVE MFD CONVERGENCE STUDY  (%s)\n', method_label);
fprintf('============================================================\n');

Nh_list  = [8, 16, 32, 64, 128, 256, 512];
n_cases  = length(Nh_list);
h_list   = 1.0 ./ Nh_list;

% Sweep of classification tolerances
tau_list = [1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4];
n_tau    = length(tau_list);

% Classification linear field: p_lin = x + y + 1
a_lin = +1;  b_lin = +1;  c_lin = 0;  d_lin = 1;

% Storage  [n_cases x n_tau]
rel_p_errors = zeros(n_cases, n_tau);
rel_m_errors = zeros(n_cases, n_tau);
tpfa_fracs   = zeros(n_cases, n_tau);
max_ind_orth = zeros(n_cases, 1);
max_ind_dist = zeros(n_cases, 1);

outDir = 'output_convergence';
if ~exist(outDir, 'dir'), mkdir(outDir); end

for it = 1:n_cases

    Nh = Nh_list(it);
    fprintf('\n============================================================\n');
    fprintf('  h = 1/%d  |  cells = %d\n', Nh, Nh^2);
    fprintf('============================================================\n');

    %% 1. Build skewed hexahedral mesh (unit cube)
    [V2, cells2D] = buildSkewed2DMesh(Nh, Nh);
    [cell_struct, face_struct, V3, ~] = extend2Dto3D(V2, cells2D, 1.0);
    n_cells = length(cell_struct);
    n_faces = length(face_struct);

    %% 2. Assign isotropic permeability K = I to all cells
    for c = 1:n_cells
        cell_struct(c).K   = eye(3);
        cell_struct(c).phi = 1.0;
        cell_struct(c).rho = 1.0;
    end

    %% 3. Build TPFA M matrix + B matrix (used in classification)
    cell_struct = createMmatrix(cell_struct, face_struct, 'tpfa');
    cell_struct = createBmatrix(cell_struct);

    face_centers     = reshape([face_struct.center], 3, [])';
    cell_centers_mat = reshape([cell_struct.center], 3, [])';

    %% ================================================================
    %% PHASE 1 — CLASSIFICATION  (GA or LA, computed ONCE per mesh)
    %% ================================================================
    d_lin_faces = a_lin*face_centers(:,1) + b_lin*face_centers(:,2) + ...
                  c_lin*face_centers(:,3) + d_lin;
    p_lin_cells = a_lin*cell_centers_mat(:,1) + b_lin*cell_centers_mat(:,2) + ...
                  c_lin*cell_centers_mat(:,3) + d_lin;

    gradp_lin = [a_lin; b_lin; c_lin];
    m_lin = zeros(n_faces, 1);
    for f = 1:n_faces
        nf     = face_struct(f).normal(:);
        m_lin(f) = -face_struct(f).area * dot(gradp_lin, nf / norm(nf));
    end

    if use_GA
        % --- Global Adaptation: accumulate residual onto faces ---
        res_3D = zeros(n_faces, 1);
        for cn = 1:n_cells
            face_ids = cell_struct(cn).faces;
            signs    = cell_struct(cn).faces_orientation(:);
            M_K  = signs .* cell_struct(cn).M;
            B_K  = cell_struct(cn).B;
            mK   = signs .* m_lin(face_ids);
            pK   = p_lin_cells(cn);
            d_K  = signs .* d_lin_faces(face_ids);
            DeltaP_K = -B_K * pK + d_K;
            R_K = (M_K * mK - B_K * pK + d_K) / norm(DeltaP_K);
            res_3D(face_ids) = res_3D(face_ids) + R_K;
        end
        % Cell indicator = max |assembled face residual| over cell faces
        cell_indicator = zeros(n_cells, 1);
        for cn = 1:n_cells
            cell_indicator(cn) = max(abs(res_3D(cell_struct(cn).faces)));
        end

    else
        % --- Local Adaptation: cell-wise inf-norm of local residual ---
        cell_indicator = zeros(n_cells, 1);
        for cn = 1:n_cells
            face_ids = cell_struct(cn).faces;
            signs    = cell_struct(cn).faces_orientation(:);
            M_K  = signs .* cell_struct(cn).M;
            B_K  = cell_struct(cn).B;
            [mK_loc, pK_loc, d_K_loc] = projectLocalAnalyticalField3D( ...
                cn, cell_struct, face_struct, face_centers, a_lin, b_lin, c_lin, d_lin);
            mK  = signs .* mK_loc;
            pK  = p_lin_cells(cn);
            d_K = signs .* d_K_loc;
            DeltaP_K = -B_K * pK + d_K;
            R_K = (M_K * mK - B_K * pK + d_K) / norm(DeltaP_K);
            cell_indicator(cn) = norm(R_K, Inf);
        end
    end

    % Distorted zone ground truth
    in_vert  = (cell_centers_mat(:,1) > 0.4 & cell_centers_mat(:,1) < 0.6);
    in_horiz = (cell_centers_mat(:,2) > 0.4 & cell_centers_mat(:,2) < 0.6);
    in_fault = in_vert | in_horiz;
    orth_vals = cell_indicator(~in_fault);
    dist_vals = cell_indicator( in_fault);
    max_ind_orth(it) = max(orth_vals);
    if isempty(dist_vals)
        max_ind_dist(it) = NaN;
        warning('No cells found in fault ribbon for Nh=%d', Nh);
    else
        max_ind_dist(it) = max(dist_vals);
    end
    fprintf('  Max indicator – K-orthogonal : %.4e\n', max_ind_orth(it));
    fprintf('  Max indicator – skewed ribbon : %.4e\n', max_ind_dist(it));

    % Exact solution (quadratic manufactured field)
    xc = cell_centers_mat(:,1);  yc = cell_centers_mat(:,2);
    p_exact = xc.^2 + 2*yc.^2 + xc.*yc;

    bnd_faces = find(arrayfun(@(f) length(f.cells) == 1, face_struct));

    m_exact = zeros(n_faces, 1);
    for f = 1:n_faces
        xf = face_struct(f).center(:);
        nf = face_struct(f).normal(:);
        u  = -[2*xf(1)+xf(2); 4*xf(2)+xf(1); 0];
        m_exact(f) = face_struct(f).area * dot(u, nf);
    end

    %% ================================================================
    %% PHASE 2 — LOOP OVER TAU  (classify -> assemble -> solve)
    %% ================================================================
    for itau = 1:n_tau
        tau = tau_list(itau);

        cellMarking = double(cell_indicator > tau);
        n_tpfa = n_cells - sum(cellMarking);
        tpfa_fracs(it, itau) = n_tpfa / n_cells;
        fprintf('  tau = %.1e | TPFA = %d / %d  (%.1f%%)\n', tau, n_tpfa, n_cells, tpfa_fracs(it,itau)*100);

        % Assemble adaptive M
        dim = 3;
        total_nnz = sum(arrayfun(@(cc) length(cc.faces)^2, cell_struct));
        rows_m = zeros(total_nnz,1);
        cols_m = zeros(total_nnz,1);
        vals_m = zeros(total_nnz,1);
        idx = 0;

        for cc = 1:n_cells
            face_ids = cell_struct(cc).faces;
            cell_nf  = length(face_ids);
            Cc       = cell_struct(cc).center(:);
            K        = cell_struct(cc).K;
            v        = cell_struct(cc).volume;
            signs    = cell_struct(cc).faces_orientation(:);
            Cf_mat   = reshape([face_struct(face_ids).center], dim, cell_nf)';
            Nf_mat   = reshape([face_struct(face_ids).normal], dim, cell_nf)';
            Af_vec   = [face_struct(face_ids).area]';
            C        = Cf_mat - Cc';
            df_norms = sqrt(sum(C.^2, 2));
            signf_vec = sign(sum((C ./ df_norms) .* Nf_mat, 2));
            N        = Af_vec .* signf_vec .* Nf_mat;
            af       = Af_vec;

            if cellMarking(cc) == 0
                td   = sum(C .* (N*K), 2) ./ sum(C .* C, 2);
                invT = diag(1 ./ abs(td));
            else
                t_loc = 6 * sum(diag(K)) / dim;
                Q     = orth(N ./ af);
                U     = eye(cell_nf) - Q*Q';
                di    = diag(1./af);
                invT  = (C*(K\C'))/v + (v/t_loc)*(di*U*di);
            end

            sign_mat = signs * signs';
            [gi, gj] = ndgrid(face_ids, face_ids);
            n2 = cell_nf^2;
            rows_m(idx+1:idx+n2) = gi(:);
            cols_m(idx+1:idx+n2) = gj(:);
            vals_m(idx+1:idx+n2) = reshape(sign_mat .* invT, [], 1);
            idx = idx + n2;
        end

        M_adapt = sparse(rows_m, cols_m, vals_m, n_faces, n_faces);
        B       = buildBmatrix(cell_struct, face_struct);

        for f = 1:n_faces
            face_struct(f).BC_pressure = [];
            face_struct(f).BC_flux     = [];
        end
        for fi = 1:length(bnd_faces)
            f  = bnd_faces(fi);
            xf = face_struct(f).center;
            face_struct(f).BC_pressure = xf(1)^2 + 2*xf(2)^2 + xf(1)*xf(2);
        end

        rhs_D = dirichletBoundary(cell_struct, face_struct);

        f_src = zeros(n_cells, 1);
        for c = 1:n_cells
            f_src(c) = 6.0 * cell_struct(c).volume;
        end

        RHS   = [rhs_D; f_src];
        A_sys = [M_adapt, -B'; B, sparse(n_cells, n_cells)];
        [~, pin_cell] = min(vecnorm(cell_centers_mat, 2, 2));
        pin_dof = n_faces + pin_cell;
        [A_bc, rhs_bc] = enforcePrescribedDOFsStrong(pin_dof, p_exact(pin_cell), A_sys, -RHS);
        sol = A_bc \ rhs_bc;

        m_num = sol(1:n_faces);
        p_num = sol(n_faces+1:end);

        rel_p_errors(it, itau) = norm(p_num - p_exact) / norm(p_exact);
        rel_m_errors(it, itau) = norm(m_num - m_exact) / norm(m_exact);
        fprintf('    Rel L2 pressure error = %.6e\n', rel_p_errors(it, itau));
        fprintf('    Rel L2 flux error     = %.6e\n', rel_m_errors(it, itau));

        vtu_data = struct(...
            'cellMarking', cellMarking, ...
            'ground_truth', double(in_fault), ...
            'indicator',   cell_indicator, ...
            'p_num',       p_num, ...
            'p_exact',     p_exact, ...
            'p_error',     abs(p_num - p_exact) ...
        );
        writeExtrudedMeshVTP( ...
            fullfile(outDir, sprintf('%s_tau%d_h%d.vtu', method_tag, itau, Nh)), ...
            V3, cell_struct, face_struct, vtu_data);
    end
end

%% ===== PLOT SETTINGS =====
fig_w = 700;  fig_h = 560;
ax_fs = 13;   title_fs = 14;  leg_fs = 10;
lw    = 2.0;  msz = 7;
markers = {'o','s','d','^','v'};
tau_colors = lines(n_tau);
h_colors   = lines(n_cases);

%% ===== PLOT 1: PRESSURE L2 CONVERGENCE =====
fig1 = figure('Position', [100 100 fig_w fig_h], 'Color', 'w');
ax1  = axes('Parent', fig1);
for itau = 1:n_tau
    loglog(ax1, h_list, rel_p_errors(:, itau), ...
        ['-' markers{mod(itau-1,numel(markers))+1}], ...
        'LineWidth', lw, 'Color', tau_colors(itau,:), 'MarkerSize', msz, ...
        'MarkerFaceColor', tau_colors(itau,:), ...
        'DisplayName', sprintf('$\\tau = 10^{%d}$', round(log10(tau_list(itau)))));
    hold(ax1, 'on');
end
ref2 = (h_list / h_list(1)).^2 * min(rel_p_errors(1,:)) * 0.15;
loglog(ax1, h_list, ref2, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
text(h_list(3), ref2(3)*0.4, '$\mathcal{O}(h^2)$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Color', [0.4 0.4 0.4], 'Parent', ax1);
set(ax1, 'FontSize', ax_fs, 'Box', 'on', 'LineWidth', 1.0, ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.25, 'MinorGridAlpha', 0.1, ...
    'TickLabelInterpreter', 'latex');
xlabel(ax1, 'Cell size $h$',           'Interpreter', 'latex', 'FontSize', ax_fs+1);
ylabel(ax1, 'Relative $L^2$ error',    'Interpreter', 'latex', 'FontSize', ax_fs+1);
title(ax1,  sprintf('Pressure Convergence  (%s)', method_label), 'FontSize', title_fs);
lgd1 = legend(ax1, 'Interpreter', 'latex', 'Location', 'southoutside', ...
    'FontSize', leg_fs, 'NumColumns', n_tau, 'Box', 'off');
exportgraphics(fig1, fullfile(outDir, sprintf('%s_pressure_convergence.pdf', method_tag)), 'ContentType', 'vector');

%% ===== PLOT 2: FLUX L2 CONVERGENCE =====
fig2 = figure('Position', [150 100 fig_w fig_h], 'Color', 'w');
ax2  = axes('Parent', fig2);
for itau = 1:n_tau
    loglog(ax2, h_list, rel_m_errors(:, itau), ...
        ['-' markers{mod(itau-1,numel(markers))+1}], ...
        'LineWidth', lw, 'Color', tau_colors(itau,:), 'MarkerSize', msz, ...
        'MarkerFaceColor', tau_colors(itau,:), ...
        'DisplayName', sprintf('$\\tau = 10^{%d}$', round(log10(tau_list(itau)))));
    hold(ax2, 'on');
end
ref1 = (h_list / h_list(1)) * min(rel_m_errors(1,:)) * 0.15;
loglog(ax2, h_list, ref1, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
text(h_list(3), ref1(3)*0.4, '$\mathcal{O}(h)$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Color', [0.4 0.4 0.4], 'Parent', ax2);
set(ax2, 'FontSize', ax_fs, 'Box', 'on', 'LineWidth', 1.0, ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.25, 'MinorGridAlpha', 0.1, ...
    'TickLabelInterpreter', 'latex');
xlabel(ax2, 'Cell size $h$',           'Interpreter', 'latex', 'FontSize', ax_fs+1);
ylabel(ax2, 'Relative $L^2$ error',    'Interpreter', 'latex', 'FontSize', ax_fs+1);
title(ax2,  sprintf('Flux Convergence  (%s)', method_label), 'FontSize', title_fs);
lgd2 = legend(ax2, 'Interpreter', 'latex', 'Location', 'southoutside', ...
    'FontSize', leg_fs, 'NumColumns', n_tau, 'Box', 'off');
exportgraphics(fig2, fullfile(outDir, sprintf('%s_flux_convergence.pdf', method_tag)), 'ContentType', 'vector');

%% ===== PLOT 3: TPFA FRACTION vs TAU =====
fig3 = figure('Position', [200 100 fig_w fig_h], 'Color', 'w');
ax3  = axes('Parent', fig3);
for it = 1:n_cases
    semilogx(ax3, tau_list, tpfa_fracs(it,:)*100, ...
        ['-' markers{mod(it-1,numel(markers))+1}], ...
        'LineWidth', lw, 'Color', h_colors(it,:), 'MarkerSize', msz, ...
        'MarkerFaceColor', h_colors(it,:), ...
        'DisplayName', sprintf('$h = 1/%d$', Nh_list(it)));
    hold(ax3, 'on');
end
set(ax3, 'FontSize', ax_fs, 'Box', 'on', 'LineWidth', 1.0, 'XDir', 'reverse', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.25, ...
    'TickLabelInterpreter', 'latex');
ylim(ax3, [0 105]);
xlabel(ax3, 'Tolerance $\tau$',        'Interpreter', 'latex', 'FontSize', ax_fs+1);
ylabel(ax3, 'TPFA fraction (\%)',      'Interpreter', 'latex', 'FontSize', ax_fs+1);
title(ax3,  sprintf('TPFA Cell Fraction vs $\\tau$  (%s)', method_label), ...
    'FontSize', title_fs, 'Interpreter', 'latex');
lgd3 = legend(ax3, 'Interpreter', 'latex', 'Location', 'southoutside', ...
    'FontSize', leg_fs, 'NumColumns', n_cases, 'Box', 'off');
exportgraphics(fig3, fullfile(outDir, sprintf('%s_tpfa_fraction.pdf', method_tag)), 'ContentType', 'vector');

%% ===== LOCAL SUBFUNCTIONS =====

function [V2, cells2D] = buildSkewed2DMesh(Nx, Ny)
% Cartesian grid on [0,1]^2 with two full-domain shear ribbons.
%   Vertical   ribbon (x in [0.4,0.6]): displaces y.  bell(x)*sin(pi*y)
%   Horizontal ribbon (y in [0.4,0.6]): displaces x.  bell(y)*sin(pi*x)
% Both shears are applied everywhere they are active (including the
% intersection square) so the displacement field is C0-continuous.
% The amplitude is kept small enough that even the worst-case compounding
% at (0.5,0.5) never produces self-intersecting quads.

    x0 = linspace(0, 1, Nx+1);
    y0 = linspace(0, 1, Ny+1);
    [X0, Y0] = meshgrid(x0, y0);
    X = X0;  Y = Y0;

    shear_amp = 0.03;   % safe for Nh >= 8: max displacement ~24% of cell size

    for i = 1:numel(X0)
        xv = X0(i);  yv = Y0(i);
        dx = 0;  dy = 0;

        % Vertical ribbon: bell in x, sine in y
        if xv >= 0.4 && xv <= 0.6
            bell_x = min(xv - 0.4, 0.6 - xv) / 0.1;   % 0->1->0
            dy = shear_amp * bell_x * sin(pi * yv);
        end

        % Horizontal ribbon: bell in y, sine in x
        if yv >= 0.4 && yv <= 0.6
            bell_y = min(yv - 0.4, 0.6 - yv) / 0.1;   % 0->1->0
            dx = shear_amp * bell_y * sin(pi * xv);
        end

        X(i) = xv + dx;
        Y(i) = yv + dy;
    end

    V2 = [X(:), Y(:)];
    cells2D = cell(Nx*Ny, 1);
    idx = 1;
    for j = 1:Ny
        for i = 1:Nx
            v1 = (j-1)*(Nx+1) + i;
            v2 = v1 + 1;
            v3 = j*(Nx+1) + i + 1;
            v4 = j*(Nx+1) + i;
            cells2D{idx} = [v1, v2, v3, v4];
            idx = idx + 1;
        end
    end
end

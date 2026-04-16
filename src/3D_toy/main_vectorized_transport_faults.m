close all; clear; clc; 
addpath(genpath('FACTORIZE'))

%% ===== Step 1: build mesh =====
setupMRSTAuto();
G = scratch_conformal();
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);
Lx = 1; Ly = 1; Lz = 0.4;

%% ===== Step 2: physical setup =====
ip_type = 'tpfa';
tol_list = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9];
n_tol = length(tol_list);

eps_solver = 1e-11;
gmres_niter = 1000;
g_c = 0.0;
dt_pressure = 1.0;

tEnd = 0.25;
dt_transport = 0.01;
Sw0 = zeros(length(cell_struct),1);
Sw_inj = 1.0;

[cell_struct, face_struct, phys] = initPhysicalParams3D(cell_struct, face_struct, Lx, Ly, Lz, 'linear');

% analytical projection
a = -1/Lx; b = 0; c = 0; d = 1;
[m_proj, p_proj] = projectAnalyticalField3D(cell_struct, face_struct, phys, a, b, c, d);

cell_struct = createMmatrix(cell_struct, face_struct, ip_type);
cell_struct = createBmatrix(cell_struct);

n_cells = length(cell_struct);
n_faces = length(face_struct);

face_centers = reshape([face_struct.center], 3, [])';
d_all = a*face_centers(:,1) + b*face_centers(:,2) + c*face_centers(:,3) + d;

%% ===== FULL MFD REFERENCE =====
dim = length(face_struct(1).center);
total_nnz = sum(arrayfun(@(cc) length(cc.faces)^2, cell_struct));

rows = zeros(total_nnz,1);
cols = zeros(total_nnz,1);
vals = zeros(total_nnz,1);
idx = 0;

for cc = 1:n_cells
    face_ids = cell_struct(cc).faces;
    cell_nf  = length(face_ids);
    Cc       = cell_struct(cc).center(:);
    K        = cell_struct(cc).K;
    v        = cell_struct(cc).volume;
    signs    = cell_struct(cc).faces_orientation(:);

    Cf_mat = reshape([face_struct(face_ids).center], dim, cell_nf)';
    Nf_mat = reshape([face_struct(face_ids).normal], dim, cell_nf)';
    Af_vec = [face_struct(face_ids).area]';

    C = Cf_mat - Cc';
    df_norms = sqrt(sum(C.^2, 2));
    signf_vec = sign(sum((C ./ df_norms) .* Nf_mat, 2));
    N = Af_vec .* signf_vec .* Nf_mat;

    % Simple
    t_loc = 6 * sum(diag(K)) / dim;
    Q  = orth(N ./ Af_vec);
    U  = eye(cell_nf) - Q * Q';
    di = diag(1 ./ Af_vec);
    invT_reg = (v / t_loc) * (di * U * di);
    invT = (C * (K \ C')) / v + invT_reg;

    % Quasi TPFA
    % W  = N * K * N';
    % Qc = orth(C);
    % P  = eye(cell_nf) - Qc * Qc';
    % diW = diag(1 ./ diag(W));
    % invT = (C * (K \ C')) / v + (v / 2) * (P * diW * P);

    % BdVLM (I NEED TO DOUBLE CHECK)
    % R  = diag(Af_vec) * C;
    % Nbd = Nf_mat * K;
    % M0 = R * ((R' * Nbd) \ R');
    % PN = eye(cell_nf) - Nbd * ((Nbd' * Nbd) \ Nbd');
    % invT = diag(1./Af_vec) * (M0 + (1/cell_nf)*PN) * diag(1./Af_vec);

    sign_mat = signs * signs';
    [gi, gj] = ndgrid(face_ids, face_ids);
    n2 = cell_nf^2;

    rows(idx+1:idx+n2) = gi(:);
    cols(idx+1:idx+n2) = gj(:);
    vals(idx+1:idx+n2) = reshape(sign_mat .* invT, [], 1);
    idx = idx + n2;
end

M_full = sparse(rows, cols, vals, n_faces, n_faces);

B = buildBmatrix(cell_struct, face_struct);
T = buildTmatrix(cell_struct);

A_full = [M_full, -B';
          B, (1/dt_pressure)*T*0];

matrix_full = A_full;

rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
[A_full, rhs_Neumann, ~] = neumannBoundary(A_full, [], face_struct);

p_n = zeros(n_cells,1);
f_g = buildGravityRHS(face_struct, g_c);

BC_face_flux.ids  = find(~arrayfun(@(s) isempty(s.BC_flux), face_struct));
BC_face_flux.vals = [face_struct(BC_face_flux.ids).BC_flux];

RHS_full = [f_g + rhs_Dirichlet;
            (1/dt_pressure)*(T*p_n)];

[matrix_full, RHS_full] = enforcePrescribedDOFsStrong( ...
    BC_face_flux.ids, BC_face_flux.vals, matrix_full, RHS_full);

sol_full = matrix_full \ (-RHS_full);
m_full = sol_full(1:n_faces);
p_full = sol_full(n_faces+1:end);

[Sw_hist_ref, time_hist] = runSinglePhaseTransportFixedFluxImplicit( ...
    cell_struct, face_struct, m_full, Sw0, Sw_inj, tEnd, dt_transport);

% Animation for saturation
% figure;
% 
% for tid = 1:length(time_hist)
%     Sw_now = Sw_hist_ref(:, tid);
% 
%     plotCellData(G, Sw_now);
%     colorbar;
%     caxis([0 1]);
% 
%     title(sprintf('t = %.3f', time_hist(tid)));
%     view(3);
% 
%     drawnow;
%     pause(0.2);
% end

Sw_ref = Sw_hist_ref(:, end);

%% ===== ADAPTIVE LOOP =====
rel_flux_errors = zeros(n_tol,1);
abs_flux_errors = zeros(n_tol,1);
sat_error_all   = zeros(n_tol,1);
abs_sat_errors  = zeros(n_tol,1);

for itol = 1:n_tol

    tol = tol_list(itol);

    res_3D = zeros(n_faces,1);

    for cn = 1:n_cells
        face_ids = cell_struct(cn).faces;
        signs    = cell_struct(cn).faces_orientation(:);

        M_K = signs .* cell_struct(cn).M;
        B_K = cell_struct(cn).B;

        mK  = signs .* m_proj(face_ids);
        pK  = p_proj(cn);
        d_K = signs .* d_all(face_ids);

        DeltaP_K = -B_K * pK + d_K;
        R_K = (M_K * mK - B_K * pK + d_K) / norm(DeltaP_K);

        res_3D(face_ids) = res_3D(face_ids) + R_K;
    end

    face_exceeds = abs(res_3D) > tol;

    face_counts  = arrayfun(@(cc) length(cc.faces), cell_struct(:));
    all_face_ids = cell2mat(arrayfun(@(cc) cc.faces(:), cell_struct(:), 'UniformOutput', false));
    all_cell_ids = repelem((1:n_cells)', face_counts);

    cellMarking_3D = accumarray(all_cell_ids, face_exceeds(all_face_ids), [n_cells 1], @any);
    tpfa_count = n_cells - sum(cellMarking_3D);
    fprintf('tol = %.1e | TPFA cells = %d / %d\n', tol, tpfa_count, n_cells);


     % Export to VTU file and save the file (only for the first tol to avoid many files)
     outDir = 'output_fault';
    
     if ~exist(outDir, 'dir')
         mkdir(outDir);
     end
    
     filename = fullfile(outDir, sprintf('mesh_l_%d.vtu', itol-1));
     writeExtrudedMeshVTP(filename, V3, cell_struct, face_struct, cellMarking_3D, 'cellMarking', 'cell_plot');

    %% pressure solve (ORIGINAL)
    rows = zeros(total_nnz,1);
    cols = zeros(total_nnz,1);
    vals = zeros(total_nnz,1);
    idx = 0;

    for cc = 1:n_cells
        face_ids = cell_struct(cc).faces;
        cell_nf  = length(face_ids);
        Cc       = cell_struct(cc).center(:);
        K        = cell_struct(cc).K;
        v        = cell_struct(cc).volume;
        signs    = cell_struct(cc).faces_orientation(:);

        Cf_mat = reshape([face_struct(face_ids).center], dim, cell_nf)';
        Nf_mat = reshape([face_struct(face_ids).normal], dim, cell_nf)';
        Af_vec = [face_struct(face_ids).area]';

        C = Cf_mat - Cc';
        df_norms = sqrt(sum(C.^2, 2));
        signf_vec = sign(sum((C ./ df_norms) .* Nf_mat, 2));
        N = Af_vec .* signf_vec .* Nf_mat;

        if cellMarking_3D(cc) == 0
            td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
            invT = diag(1 ./ abs(td));
        else
            % Simple
            t_loc = 6 * sum(diag(K)) / dim;
            Q  = orth(N ./ Af_vec);
            U  = eye(cell_nf) - Q * Q';
            di = diag(1 ./ Af_vec);
            invT_reg = (v / t_loc) * (di * U * di);
            invT = (C * (K \ C')) / v + invT_reg;

            % Quasi TPFA
            % W  = N * K * N';
            % Qc = orth(C);
            % P  = eye(cell_nf) - Qc * Qc';
            % diW = diag(1 ./ diag(W));
            % invT = (C * (K \ C')) / v + (v / 2) * (P * diW * P);

            % BdVLM (Again, I need to double check!)
            % R  = diag(Af_vec) * C;
            % Nbd = Nf_mat * K;
            % M0 = R * ((R' * Nbd) \ R');
            % PN = eye(cell_nf) - Nbd * ((Nbd' * Nbd) \ Nbd');
            % invT = diag(1./Af_vec) * (M0 + (1/cell_nf)*PN) * diag(1./Af_vec);
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

    if itol == 1
        nnz_full = nnz(M_full);
    end

    nnz_adaptive = nnz(M);

    sparsity_reduction = 1 - nnz_adaptive / nnz_full;

    fprintf('nnz(full MFD) = %d, nnz(adaptive) = %d, reduction = %.2f%%\n', ...
        nnz_full, nnz_adaptive, 100*sparsity_reduction);

    A_full = [M, -B';
              B, (1/dt_pressure)*T*0];

    if abs(tol - 1e-5) < 1e-12
        A_tau = A_full;
    end

    matrix = A_full;

    [matrix, RHS] = enforcePrescribedDOFsStrong( ...
        BC_face_flux.ids, BC_face_flux.vals, matrix, RHS_full);

    sol3 = matrix \ (-RHS);

    m_num = sol3(1:n_faces);
    rel_flux_errors(itol) = norm(m_num - m_full) / norm(m_full);
    abs_flux_errors(itol) = norm(m_num - m_full);

    plot_flux_error = max(rel_flux_errors, 1e-16);
    plot_abs_flux_error = max(abs_flux_errors, 1e-16);

    %% transport
    [Sw_hist, ~] = runSinglePhaseTransportFixedFluxImplicit( ...
        cell_struct, face_struct, m_num, Sw0, Sw_inj, tEnd, dt_transport);

    Sw_final = Sw_hist(:, end);
    Sw_all{itol} = Sw_hist;
    sat_error_all(itol)  = norm(Sw_final - Sw_ref) / norm(Sw_ref);
    abs_sat_errors(itol) = norm(Sw_final - Sw_ref);

    plot_sat_error = max(sat_error_all, 1e-16);
    plot_abs_sat_error = max(abs_sat_errors, 1e-16);

    fprintf('tol = %.1e | sat error = %.6e\n', tol, sat_error_all(itol));
    fprintf('tol = %.1e | sat error = %.6e\n', tol, abs_sat_errors(itol));

    filename = fullfile(outDir, sprintf('sat_tol_%d.vtu', itol));
    writeExtrudedMeshVTP(filename, V3, cell_struct, face_struct, Sw_final, 'saturation', 'saturation_plot');
end

Sw1 = Sw_all{1}(:, end); 
Sw2 = Sw_all{10}(:, end);

%% ===== FLUX ERROR PLOT =====
figure;

c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];

loglog(tol_list, plot_flux_error, '-o', 'LineWidth', 2, 'Color', c1, ...
    'DisplayName', 'Relative flux error');
hold on;

loglog(tol_list, plot_abs_flux_error, '-s', 'LineWidth', 2, 'Color', c2, ...
    'DisplayName', 'Absolute flux error');

loglog(tol_list, tol_list, '--', 'LineWidth', 1.5, 'Color', c1, ...
    'DisplayName', 'y = tol');

loglog(tol_list, tol_list, ':', 'LineWidth', 1.5, 'Color', c2, ...
    'DisplayName', 'y = tol');

xlabel('Tolerance', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
title('Mass Flux Error vs Tolerance', 'FontSize', 18);

legend('Location', 'northwest', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([min(tol_list) max(tol_list)]);
grid on;

set(gcf, 'Units','pixels','Position',[100 100 1200 900]);

%% ===== SATURATION PLOT =====
figure;

loglog(tol_list, plot_sat_error, '-o', 'LineWidth', 2, 'Color', c1, ...
    'DisplayName', 'Relative saturation error');
hold on;

loglog(tol_list, plot_abs_sat_error, '-s', 'LineWidth', 2, 'Color', c2, ...
    'DisplayName', 'Absolute saturation error');

loglog(tol_list, tol_list, '--', 'LineWidth', 1.5, 'Color', c1, ...
    'DisplayName', 'y = tol');

loglog(tol_list, tol_list, ':', 'LineWidth', 1.5, 'Color', c2, ...
    'DisplayName', 'y = tol');

xlabel('Tolerance', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
title('Saturation Error vs Tolerance', 'FontSize', 18);

legend('Location', 'northwest', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([min(tol_list) max(tol_list)]);
grid on;

set(gcf, 'Units','pixels','Position',[100 100 1200 900]);

%%%%%%%%%%%%%%%%% HELPER FUNCTION!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sw_hist, time_hist] = runSinglePhaseTransportFixedFluxImplicit( ...
    cell_struct, face_struct, m_num, Sw0, Sw_inj, tEnd, dt)

    n_cells = length(cell_struct);
    
    Sw = Sw0(:);
    t = 0;
    
    Sw_hist = Sw;
    time_hist = t;
    
    Vc  = arrayfun(@(c) c.volume, cell_struct)';
    phi = arrayfun(@(c) c.phi, cell_struct)';
    acc = phi .* Vc;
    
    while t < tEnd
        dt_step = min(dt, tEnd - t);
    
        rows = [];
        cols = [];
        vals = [];
        rhs  = (acc / dt_step) .* Sw;
    
        for c = 1:n_cells
            rows(end+1,1) = c;
            cols(end+1,1) = c;
            vals(end+1,1) = acc(c) / dt_step;
    
            faces = cell_struct(c).faces;
            sgns  = cell_struct(c).faces_orientation;
    
            for k = 1:length(faces)
                f = faces(k);
                Fcf = sgns(k) * m_num(f);
                neigh = face_struct(f).cells;
    
                if numel(neigh) == 2
                    other = neigh(neigh ~= c);
    
                    if Fcf >= 0
                        rows(end+1,1) = c;
                        cols(end+1,1) = c;
                        vals(end+1,1) = Fcf;
                    else
                        rows(end+1,1) = c;
                        cols(end+1,1) = other;
                        vals(end+1,1) = Fcf;
                    end
    
                elseif numel(neigh) == 1
                    if Fcf >= 0
                        rows(end+1,1) = c;
                        cols(end+1,1) = c;
                        vals(end+1,1) = Fcf;
                    else
                        rhs(c) = rhs(c) - Fcf * Sw_inj;
                    end
                end
            end
        end
    
        A = sparse(rows, cols, vals, n_cells, n_cells);
    
        Sw = A \ rhs;
        Sw = max(0, min(1, Sw));
    
        t = t + dt_step;
        Sw_hist(:, end+1) = Sw;
        time_hist(end+1) = t;
    end
end
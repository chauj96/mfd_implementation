% clear all; clc; 
addpath('PolyMesher/');

% --- Basic settings ---
case_type = 'structured';
dt = 1;
nx = 20; nz = 20;
Lx = 1.0; Lz = 1.0;
rho = 1000;
g_c = 0.0;
ip_type = 'tpfa';

% --- tolerance sweep values ---
% tols = 1e-16;
tols = 1e-4;
% tols = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15];


% --- storage arrays ---
flux_error = zeros(size(tols));
tpfa_ratio = zeros(size(tols));
mean_logR  = zeros(size(tols));

% --- Build grid once (no distortion yet) ---
if strcmp(case_type, 'structured')
    [cell_struct_base, face_struct_base, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz);
else
    domain = @MbbDomain;
    n_cells = 100;
    [cell_struct_base, face_struct_base, vertices, cells] = buildPolyGrid(domain, n_cells);
end

dummy = 0;
% Sweep over tolerance values
for k = 1:length(tols)
    tol_R = tols(k);
    fprintf('\n=== Running case tol_R = %.1e ===\n', tol_R);

    % copy base grid
    cell_struct = cell_struct_base;
    face_struct = face_struct_base;

    % Physical initialization
    [cell_struct, face_struct, phys] = initPhysicalParams(cell_struct, face_struct, Lx, Lz, case_type);

    % Build matrices
    [M, Rvals, scheme_labels, cell_struct] = buildMmatrixParametric(cell_struct, face_struct, ip_type, tol_R);
    B = buildBmatrix(cell_struct, face_struct);
    T = buildTmatrix(cell_struct);
    A = [M, -B'; B, (1/dt)*T];

    % Boundary conditions
    rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
    [A, rhs_Neumann] = neumannBoundary(A, cell_struct, face_struct);

    % RHS
    p_n = zeros(length(cell_struct),1);
    f_g = buildGravityRHS(face_struct, g_c);
    rhs = [f_g; (1/dt)*(T*p_n)];
    rhs_BC = [rhs_Dirichlet + rhs_Neumann; zeros(length(cell_struct),1)];

    sol = A \ -(rhs + rhs_BC);
    n_faces = length(face_struct);
    m_sol = sol(1:n_faces);
    p_sol = sol(n_faces+1:end);

    % Flux L2 error
    flux_error(k) = sqrt(sum((m_sol - phys.m_ref_faces).^2) / sum(phys.m_ref_faces.^2));
   
    nTPFA = sum(strcmp(scheme_labels,'T'));
    tpfa_ratio(k) = nTPFA;
    Rvals = abs(log(Rvals));
    mean_logR(k) = mean(Rvals);
    dummy = Rvals;

    fprintf('Flux L2 error     : %.3e\n', flux_error(k));
    fprintf('TPFA ratio        : %.3f\n', tpfa_ratio(k));
    fprintf('Mean log10(R)     : %.12f\n', mean_logR(k));
end


% Plot only L2 error 
% hold on;
% figure;
% loglog(tols, flux_error, '-o','LineWidth',1.5,'MarkerSize',6);
% xlabel('Tolerance \it{tol_R}','FontSize',12);
% ylabel('Flux L2 relative error','FontSize',12);
% title('Flux Error vs Tolerance','FontSize',13);
% grid on;
% 
plotPressurePolygonal(vertices, cells, p_sol, ip_type, scheme_labels);
% 
% % Display summary table
fprintf('\n=== Summary Table ===\n');
fprintf('%12s %15s %15s %15s\n','tol_R','Flux L2 error','TPFA ratio','Mean log10(R)');
for k = 1:length(tols)
    fprintf('%12.1e %15.3e %15.3f %15.3f\n', ...
        tols(k), flux_error(k), tpfa_ratio(k), mean_logR(k));
end


% ====================================================================

% clear all; clc; close all
% 
% addpath('PolyMesher/');
% 
% % case options: 'structured' or 'unstructured'
% case_type = 'structured';
% 
% dt = 1;
% nx = 50;
% nz = 50;
% Lx = 1.0;
% Lz = 1.0;
% rho = 1000;
% g_c = 0.0 * (1.0e-6);
% % g_c = 0;
% 
% if strcmp(case_type, 'structured')
%     [cell_struct, face_struct, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz);
% 
% elseif strcmp(case_type, 'unstructured')
%     domain = @MbbDomain; % this domain also sets with Lx = 1, Lz = 1 -> BdBox = [0 1 0 1];
%     % n_cells = nz*nx; % we can change the number of cells
%     n_cells = 25;
%     [cell_struct, face_struct, vertices, cells] = buildPolyGrid(domain, n_cells);
% end
% 
% % Step 1: Initialize physical properties
% [cell_struct, face_struct, phys] = initPhysicalParams(cell_struct, face_struct, Lx, Lz, case_type);
% 
% % Step 2: Build matrix M, B, T / Assemble the matrices
% % TPFA case
% ip_type = 'tpfa';
% tol_R = 1e-3;
% [M, R, scheme_labels, cell_struct] = buildMmatrixParametric(cell_struct, face_struct, ip_type, tol_R);
% 
% B = buildBmatrix(cell_struct, face_struct);
% T = buildTmatrix(cell_struct);
% 
% A = [M, -B';B, (1/dt)*T];
% 
% % figure;
% % spy(A);                 % Display as colored image
% % axis equal tight;             % Equal aspect ratio, no padding
% % colormap(parula);             % Try 'parula', 'viridis', or 'turbo'
% % colorbar;                     % Add colorbar
% % % set(gca, 'YDir', 'normal');   % Make y-axis go bottom to top
% % box on;
% 
% % Step 3: Apply boundary conditions (Dirichlet, Neumann)
% rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
% [A, rhs_Neumann] = neumannBoundary(A, cell_struct, face_struct);
% 
% % Step 4: Build RHS
% % p_n = ones(length(cell_struct),1) * 1e5; % constant initial pressure
% z_top = max(arrayfun(@(c) c.center(2), cell_struct)); 
% 
% % Hydrostatic pressure profile
% p_n = zeros(length(cell_struct),1);
% % for i = 1:length(cell_struct)
% %     z_i = cell_struct(i).center(2);
% %     p_n(i) = 1e5 + rho * g * (z_top - z_i);
% % end
% 
% f_g = buildGravityRHS(face_struct, g_c);
% rhs = [f_g; (1/dt) * (T * p_n)];
% 
% rhs_BC = [rhs_Dirichlet + rhs_Neumann ; zeros(length(cell_struct),1)];
% 
% 
% % Step 5: Solve linear system
% sol = A \ -(rhs + rhs_BC);
% 
% n_faces = length(face_struct);
% m_sol = sol(1:n_faces);
% p_sol = sol(n_faces+1:end);
% 
% % flux error!!!!!
% flux_error = sqrt(sum((m_sol - phys.m_ref_faces).^2) / sum(phys.m_ref_faces.^2));
% fprintf('Flux L2 relative error: %.3e\n', flux_error);
% 
% % Step 6: Plot the pressure field
% plotPressurePolygonal(vertices, cells, p_sol, ip_type, scheme_labels);








% Step 7: Plot log(R) at cell centers as 3D markers above discretized domain
% figure; hold on;
% 
% for i = 1:numel(cells)
%     poly = vertices(cells{i}, :);
%     fill3(poly(:,1), poly(:,2), zeros(size(poly,1),1), ...
%           [0.9 0.9 0.9], 'EdgeColor', [0.6 0.6 0.6]);
% end
% 
% x = arrayfun(@(c) c.center(1), cell_struct);
% y = arrayfun(@(c) c.center(2), cell_struct);
% z = log10(arrayfun(@(c) c.R, cell_struct));
% 
% scatter3(x, y, z, 80, z, 'filled', 'MarkerEdgeColor','k');
% colormap jet; colorbar;
% xlabel('x'); ylabel('y'); zlabel('log_{10}(R)');
% title('log(R) values at cell centers');
% view(30,30);
% grid on;
% 
% for i = 1:numel(x)
%     plot3([x(i), x(i)], [y(i), y(i)], [0, z(i)], 'k-', 'LineWidth', 0.8);
% end

% Step 7: Plot stream lines from flux
% plotStreamlinesFromFlux(cell_struct, face_struct, m_sol);


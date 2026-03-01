% clear all; clc;

% You may select one of the predefined options below.
% Any custom 2D mesh can also be used, provided that it defines:
% vertices, cell-to-vertex connectivity and follows the face-based
% structure required by the solver. The subsequent 3D extrusion and
% discretization steps are independent of how the 2D mesh is generated.

%% Step 1: build 2D meshes 

% (Option 1): SPE11B Model (read CSV geometry files)
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
nx = 30;
ny = 30;
[cell2D, face2D, V2, cells2D] = buildTestMesh(nx, ny, Lx, Ly);

%% Step 2: Extrude 2D mesh along vertical direction (H = height) to get 3D mesh
[cell_struct, face_struct, V3, cells3D] = extend2Dto3D(V2, cells2D, H);

%% Step 3: Assign physical values and get projection of analytical solution (flux and pressure)
ip_type = 'tpfa';
tol = 1e-6; % Tolerance for residual-based cell classification
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

res_3D = zeros(length(face_struct),1);

for cn = 1:length(cell_struct)
    face_ids = cell_struct(cn).faces;

    M_K = cell_struct(cn).M;          
    B_K = cell_struct(cn).B;          
    d_K = computeLocalDirichlet3D(cn, cell_struct, face_struct, a, b, c, d);

    mK = m_proj(face_ids);          
    pK = p_proj(cn);                 

    R_K = M_K * mK - B_K * pK + d_K; % residual equation
    res_3D(face_ids) = res_3D(face_ids) + R_K; % mapping to global residual vector
end

cellMarking_3D = zeros(length(cell_struct),1);
for i = 1:length(cell_struct)
    arr = abs(res_3D(cell_struct(i).faces)) > tol;
    if sum(arr) >= 1
        cellMarking_3D(i,1) = 1;
    end
end
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

    % build local geometry
    C = zeros(cell_nf, dim);
    N = zeros(cell_nf, dim);
    a = zeros(cell_nf, 1);

    for k = 1:cell_nf
        f = face_ids(k);
        Cf = face_struct(f).center(:);
        Nf = face_struct(f).normal(:);
        Af = face_struct(f).area;

        df = Cf - Cc;
        signf = sign(df'/norm(df) * Nf);
        assert(signf == signs(k), 'Orientation mismatch in cell %d', c);

        C(k,:) = df';
        N(k,:) = Af * signf * Nf';
        a(k) = Af;
    end

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


    % Global assembly
    for i = 1:cell_nf
        fi = face_ids(i); si = signs(i);
        for j = 1:cell_nf
            fj = face_ids(j); sj = signs(j);
            idx = idx + 1;
            rows(idx) = fi;
            cols(idx) = fj;
            vals(idx) = si * sj * invT(i,j);
        end
    end
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

% Solve the linear system with back slash - Schur complement didn't work :(
sol = A_full \ (-rhs_total);

m_num = sol(1:n_faces);
p_num = sol(n_faces+1:end);

% Maximum error / Relative error
[maxErr, maxIdx] = max(abs(m_num - m_proj));
relErr = norm(m_num - m_proj) / norm(m_proj);
fprintf('Relative error = %e\n', relErr);

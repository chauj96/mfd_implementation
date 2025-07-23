clear all; clc; close all;


ip_types = {'tpfa', 'simple', 'general_parametric'};

epsilon_vals = linspace(0, 3/7, 5);  % test multiple distortions
tol = 1e-10;  % tolerance for exactness checks

for e_idx = 1:length(epsilon_vals)
    epsilon = epsilon_vals(e_idx);
    fprintf('\nTesting with epsilon = %.5f\n', epsilon);
    
    % Build 2-cell mesh with distortion epsilon
    [cell_struct, face_struct, vertices, cells] = buildStructured2CellMesh(epsilon);
    
    % Initialize physical properties if needed, or assign permeability to identity
    [cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct, 1.0, 1.0, 'structured');
    
    % Build fixed matrices
    B = buildBmatrix(cell_struct, face_struct);
    
    n_faces = length(face_struct);
    n_cells = length(cell_struct);

    rhs_partial = [zeros(n_faces, 1); zeros(n_cells, 1)];
    rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct);
    
    
    p_solutions = struct();
    
    for i = 1:length(ip_types)
        ip_type = ip_types{i};
        
        M = buildMmatrixParametric(cell_struct, face_struct, ip_type);
        
        A = [M, -B'; B, 0.0 * eye(n_cells)];
        [A, rhs_Neumann] = neumannBoundary(A, cell_struct, face_struct);
        rhs_BC = [rhs_Dirichlet + rhs_Neumann; zeros(length(cell_struct), 1)];
        rhs = rhs_partial + rhs_BC;
        
        sol = A \ -rhs;
        m_sol = sol(1:n_faces)
        p_sol = sol(n_faces+1:end);
        p_solutions.(ip_type) = p_sol;
        
        centers_matrix = cell2mat({cell_struct.center}'); % 2 x n_cells
        x_coords = centers_matrix(:,1);
        [m_e, p_e] = projectExactSolution(cell_struct, face_struct);
        % For epsilon=0, TPFA should be exact, check that
        if strcmp(ip_type, 'tpfa') && abs(epsilon) < 1e-14
            err_m = norm(m_sol - m_e);
            err_p = norm(p_sol - p_e);
            assert(err_m < tol, 'TPFA m solution is not exact for zero distortion!');
            assert(err_p < tol, 'TPFA p solution is not exact for zero distortion!');
        end
        
        % For simple and general_parametric, solution should be exact for all epsilon in [0, 3/7]
        if (strcmp(ip_type, 'simple') || strcmp(ip_type, 'general_parametric'))
            err_m = norm(m_sol - m_e);
            err_p = norm(p_sol - p_e);
            assert(err_p < tol && err_m < tol, ...
                sprintf('%s inner product failed consistency test at epsilon=%.5f', ip_type, epsilon));
        end
    end
    
    % Optional: display or plot
    fprintf('All IPs passed consistency for epsilon=%.5f\n', epsilon);
end


% ====== Local function for mesh building ======
function [cell_struct, face_struct, vertices, cells] = buildStructured2CellMesh(epsilon)
    if nargin < 1
        epsilon = 0;
    end

    vertices = [...
        0.0, 0.0;   % v1
        0.5 - epsilon, 0.0;   % v2
        0.5 + epsilon, 1.0; % v3        
        0.0, 1.0;   % v4
        1.0, 0.0;  % v5
        1.0, 1.0];   % v6
        

    cells = { [1, 2, 3, 4], [2, 5, 6, 3] };

    face_map = containers.Map;
    face_struct = struct('cells', {}, 'normal', {}, 'center', {}, 'area', {});
    cell_struct = struct('center', {}, 'faces', {}, 'volume', {});

    face_counter = 0;

    for c = 1:2
        verts = cells{c};
        n_verts = length(verts);

        pts = vertices(verts, :);
        cell_center = mean(pts, 1);
        x = pts(:,1); z = pts(:,2);
        area = 0.5 * abs(sum(x .* circshift(z,-1)) - sum(z .* circshift(x,-1)));

        cell_struct(c).center = cell_center;
        cell_struct(c).volume = area;

        face_ids = zeros(1, n_verts);
        faces_orientation = zeros(1,n_verts);

        v1 = verts(1);
        v2 = verts(2);
        v3 = verts(3);
        v4 = verts(4);
        edge_list = {[v1, v2], [v2, v3], [v3, v4], [v4, v1]};
        for i = 1:n_verts
            edge = edge_list{i};
            % Create an unordered key to ensure uniqueness of the face
            edge_sorted = sort(edge);
            key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));
            if isKey(face_map, key)
                f_id = face_map(key);
                face_struct(f_id).cells(end+1) = c;
            else
                face_counter = face_counter + 1;
                f_id = face_counter;
                face_map(key) = f_id;

                p1 = vertices(edge_sorted(1),:);
                p2 = vertices(edge_sorted(2),:);
                t = p2 - p1;
                normal = [-t(2), t(1)] / norm(t);
                center = 0.5 * (p1 + p2);
                length_face = norm(t);

                face_struct(f_id).cells = c;
                face_struct(f_id).normal = normal;
                face_struct(f_id).center = center;
                face_struct(f_id).area = length_face;
            end

            face_ids(i) = f_id;

            % Compute face orientation sign for this cell:
            nf = face_struct(f_id).normal;
            xf = face_struct(f_id).center(:);
            faces_orientation(i) = sign(dot(nf, xf - cell_center'));
        end

        cell_struct(c).faces = face_ids;
        cell_struct(c).faces_orientation = faces_orientation;
    end
end

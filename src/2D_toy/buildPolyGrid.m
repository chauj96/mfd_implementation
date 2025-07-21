function [cell_struct, face_struct, vertices, cells] = buildPolyGrid(domain, n_cells)
 
    % number of Lloyd iteration
    maxIter = 30; 
    [vertices, cells, ~, ~, cell_centers] = PolyMesher(domain, n_cells, maxIter);
    
    % FACE STORAGE: center, faces, face_dirs, volume
    face_map = containers.Map();
    face_count = 0;
    face_struct = struct('nodes', {}, 'cells', {}, 'center', {}, 'normal', {}, 'area', {});
    cell_struct = struct('center', {}, 'faces', {}, 'face_dirs', {}, 'volume', {});
    
    for c = 1:n_cells
        verts = cells{c}; % vertex indices (global) of cell c
        n_verts = length(verts);
        
        for i = 1:n_verts
            v1 = verts(i);
            v2 = verts(mod(i, n_verts) + 1); % wrap around
    
            % make unique not only the key but the normal computation 
            % if v2 < v1
            %     tmp = v1;
            %     v1 = v2;
            %     v2 = tmp;
            % end
            % key1 = sprintf('%d_%d', v1, v2);
            key1 = sprintf('%d_%d', min(v1, v2), max(v1, v2)); % unordered key (edge is undirected)
            if isKey(face_map, key1)
                f_id = face_map(key1);
                face_struct(f_id).cells(end+1) = c;
            else
                face_count = face_count + 1;
                face_map(key1) = face_count;
    
                % build new face
                x1 = vertices(v1,:);
                x2 = vertices(v2,:);
                center = 0.5 * (x1 + x2);
                normal = [x2(2) - x1(2); -(x2(1) - x1(1))]; % 90 degrees rotation [cos90 sin90; -sin90 cos90]
                normal = normal / norm(normal); 
                area = norm(x2 - x1);
    
                face_struct(face_count).nodes = [v1, v2];
                face_struct(face_count).cells = [c];
                face_struct(face_count).center = center(:);
                face_struct(face_count).normal = normal(:);
                face_struct(face_count).area = area;
            end
        end
    end
    
    % CELL STORAGE: center, faces, face_dirs, volume
    for c = 1:n_cells
        verts = cells{c};
        n_verts = length(verts);
    
        cell_struct(c).center = cell_centers(c,:)';
        cell_struct(c).volume = polyarea(vertices(verts,1), vertices(verts,2));
    
        cell_struct(c).faces = zeros(1, n_verts);
        cell_struct(c).face_dirs = zeros(1, n_verts);
    
        xc = cell_struct(c).center;
    
        for i = 1:n_verts
            v1 = verts(i);
            v2 = verts(mod(i, n_verts) + 1);
    
            key1 = sprintf('%d_%d', min(v1, v2), max(v1, v2));
            f_id = face_map(key1);
    
            nf = face_struct(f_id).normal;
            xf = face_struct(f_id).center;
    
            dir = sign(dot(xf - xc, nf));
    
            cell_struct(c).faces(i) = f_id;
            cell_struct(c).face_dirs(i) = dir;
        end
        
    end
end



function [cell_struct, face_struct, vertices, cells] = buildPolyGrid(domain, n_cells)
 
    % number of Lloyd iteration
    maxIter = 30; 
    [vertices, cells, ~, ~, cell_centers] = PolyMesher(domain, n_cells, maxIter);
    
    % Initialize face and cell storage
    face_map = containers.Map();
    face_count = 0;
    face_struct = struct('nodes', {}, 'cells', {}, 'center', {}, 'normal', {}, 'area', {});
    cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'volume', {});

    for c = 1:n_cells
        verts = cells{c};
        n_verts = length(verts);
        xc = cell_centers(c,:)';

        for i = 1:n_verts
            v1 = verts(i);
            v2 = verts(mod(i, n_verts) + 1); % wrap around

            edge = [v1, v2];
            edge_sorted = sort(edge);
            key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));

            if isKey(face_map, key)
                f = face_map(key);
                face_struct(f).cells(end+1) = c;
            else
                face_count = face_count + 1;
                f = face_count;
                face_map(key) = f;

                p1 = vertices(edge_sorted(1),:);
                p2 = vertices(edge_sorted(2),:);
                t = p2 - p1;

                % Compute normal by 90-degree CCW rotation
                normal = [-t(2); t(1)];
                normal = normal / norm(normal);

                center = 0.5 * (p1 + p2);
                area = norm(t);

                face_struct(f).nodes = [v1, v2];
                face_struct(f).cells = [c];
                face_struct(f).center = center(:);
                face_struct(f).normal = normal(:);
                face_struct(f).area = area;
            end
        end
    end

    % Assemble cell data
    for c = 1:n_cells
        verts = cells{c};
        n_verts = length(verts);

        xc = cell_centers(c,:)';
        cell_struct(c).center = xc;
        cell_struct(c).volume = polyarea(vertices(verts,1), vertices(verts,2));

        face_ids = zeros(1, n_verts);
        faces_orientation = zeros(1, n_verts);

        for i = 1:n_verts
            v1 = verts(i);
            v2 = verts(mod(i, n_verts) + 1);
            edge_sorted = sort([v1, v2]);
            key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));
            f_id = face_map(key);

            nf = face_struct(f_id).normal;
            xf = face_struct(f_id).center(:);
            faces_orientation(i) = sign(dot(nf, xf - xc));

            face_ids(i) = f_id;
        end

        cell_struct(c).faces = face_ids;
        cell_struct(c).faces_orientation = faces_orientation;
    end
end

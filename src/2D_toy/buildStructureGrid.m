function [cell_struct, face_struct, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz)
 
    % create vertex list
    xv = linspace(0, Lx, nx+1);
    zv = linspace(0, Lz, nz+1);
    [X, Z] = meshgrid(xv, zv);
    X = X';
    Z = Z';

    lambda_x = 5;
    lambda_z = 5;
    amp_x = 0.0 * 0.00625;
    amp_z = 0.0 * 0.00625;

    for i = 2:nx
        for j = 2:nz
            X(i,j) = X(i,j) + amp_x * sin(lambda_x*pi * X(i,j)) * sin(lambda_z*pi * Z(i,j));
            Z(i,j) = Z(i,j) + amp_z * cos(lambda_x*pi * X(i,j)) * sin(lambda_z*pi * Z(i,j));
        end
    end

    vertices = [X(:), Z(:)];

    % FACE STORAGE
    face_struct = struct('cells', {}, 'normal', {}, 'center', {}, 'area', {});
    face_counter = 0;
    face_map = containers.Map;

    % CELL STORAGE
    cell_struct = struct('center', {}, 'faces', {}, 'face_dirs', {}, 'volume', {});
    cells = {};
    cell_id = 0;

    for j = 1:nz
        for i = 1:nx
            cell_id = cell_id + 1;
            % vertex IDs (counter clockwise)
            v1 = (j-1)*(nx+1) + i; % bottom left
            v2 = v1 + 1; % bottom right
            v3 = v2 + (nx+1); % top right
            v4 = v1 + (nx+1); % top left
            vids = [v1, v2, v3, v4];
            cells{cell_id} = vids;

            % cell center and area from vertices
            pts = vertices(vids, :);
            cell_struct(cell_id).center = mean(pts, 1);
            x = pts(:,1);
            z = pts(:,2);
            area = 0.5 * abs(sum(x .* circshift(z,-1)) - sum(z .* circshift(x,-1))); % shoelace formula
            cell_struct(cell_id).volume = area;

            % add faces (edges of the polygon)
            face_ids = [];
            face_dirs = [];

            % Define cell edges (v1-v2, v2-v3, v3-v4, v4-v1)
            edge_list = {[v1, v2], [v2, v3], [v3, v4], [v4, v1]};
            for k = 1:4
                edge = sort(edge_list{k});
                key = sprintf('%d_%d', edge(1), edge(2));
                if isKey(face_map, key)
                    f = face_map(key);
                    face_struct(f).cells = [face_struct(f).cells, cell_id];
                    face_ids(end+1) = f;
                    face_dirs(end+1) = -1;
                else
                    face_counter = face_counter + 1;
                    f = face_counter;
                    face_map(key) = f;

                    p1 = vertices(edge(1),:);
                    p2 = vertices(edge(2),:);
                    t = p2 - p1;
                    normal = [-t(2); t(1)] / norm(t);
                    center = 0.5 * (p1 + p2);
                    length_face = norm(t);

                    face_struct(f).cells = cell_id;
                    face_struct(f).normal = normal;
                    face_struct(f).center = center;
                    face_struct(f).area = length_face;

                    face_ids(end+1) = f;
                    face_dirs(end+1) = 1;
                end
            end

            cell_struct(cell_id).faces = face_ids;
            cell_struct(cell_id).face_dirs = face_dirs;
        end
    end
end

function [cell_struct, face_struct, vertices, cells] = buildStructureGrid(nx, nz, Lx, Lz)
 
    % Create vertex list
    xv = linspace(0, Lx, nx+1);
    zv = linspace(0, Lz, nz+1);
    [X, Z] = meshgrid(xv, zv);
    X = X';
    Z = Z';

    % Perturb internal vertices only (exclude boundaries)
    lambda_x = 4;
    lambda_z = 4;
    amp_x = 4.0*0.0125;
    amp_z = 4.0*0.0125;
    for i = 2:nx
        for j = 2:nz
            X(i,j) = X(i,j) + amp_x * sin(lambda_x*pi * X(i,j)) * sin(lambda_z*pi * Z(i,j));
            Z(i,j) = Z(i,j) + amp_z * sin(lambda_x*pi * X(i,j)) * sin(lambda_z*pi * Z(i,j));
        end
    end

    vertices = [X(:), Z(:)];

    % Initialize face storage
    face_struct = struct('cells', {}, 'normal', {}, 'center', {}, 'area', {});
    face_counter = 0;
    face_map = containers.Map;

    % Initialize cell storage
    cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'volume', {});
    cells = {};
    cell_id = 0;

    for j = 1:nz
        for i = 1:nx
            cell_id = cell_id + 1;

            % Vertex IDs for the current cell (counter clockwise)
            v1 = (j-1)*(nx+1) + i;       % bottom left
            v2 = v1 + 1;                 % bottom right
            v3 = v2 + (nx+1);            % top right
            v4 = v1 + (nx+1);            % top left
            vids = [v1, v2, v3, v4];
            cells{cell_id} = vids;

            % Compute cell center and area using shoelace formula
            pts = vertices(vids, :);
            xc = mean(pts, 1)';
            cell_struct(cell_id).center = xc;
            x = pts(:,1);
            z = pts(:,2);
            area = 0.5 * abs(sum(x .* circshift(z,-1)) - sum(z .* circshift(x,-1)));
            cell_struct(cell_id).volume = area;

            % Process cell faces (edges)
            face_ids = zeros(1,4);
            faces_orientation = zeros(1,4);

            % Cell edges (v1-v2, v2-v3, v3-v4, v4-v1)
            edge_list = {[v1, v2], [v2, v3], [v3, v4], [v4, v1]};

            for k = 1:4
                edge = edge_list{k};
                % Create an unordered key to ensure uniqueness of the face
                edge_sorted = sort(edge);
                key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));

                if isKey(face_map, key)
                    f = face_map(key);
                    % Append this cell to face's cell list
                    face_struct(f).cells(end+1) = cell_id;
                else
                    % New face
                    face_counter = face_counter + 1;
                    f = face_counter;
                    face_map(key) = f;

                    p1 = vertices(edge_sorted(1),:);
                    p2 = vertices(edge_sorted(2),:);
                    t = p2 - p1;

                    % Compute outward normal (90 degree CCW rotation)
                    normal = [-t(2); t(1)];
                    normal = normal / norm(normal);

                    center = 0.5 * (p1 + p2);
                    length_face = norm(t);

                    face_struct(f).cells = cell_id;
                    face_struct(f).normal = normal;
                    face_struct(f).center = center;
                    face_struct(f).area = length_face;
                end

                face_ids(k) = f;

                % Compute face orientation sign for this cell:
                nf = face_struct(f).normal;
                xf = face_struct(f).center(:);
                faces_orientation(k) = sign(dot(nf, xf - xc));
            end

            cell_struct(cell_id).faces = face_ids;
            cell_struct(cell_id).faces_orientation = faces_orientation;
        end
    end
end

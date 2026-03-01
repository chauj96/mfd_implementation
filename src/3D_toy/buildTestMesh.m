function [cell_struct, face_struct, vertices, cells] = buildTestMesh(nx, nz, Lx, Lz)
% Construct a 2D structured quadrilateral mesh with a localized smooth
% geometric distortion applied inside a rotated elliptical region

    % Create structured vertices
    xv = linspace(0, Lx, nx+1);
    zv = linspace(0, Lz, nz+1);
    [X, Z] = meshgrid(xv, zv);
    X = X';
    Z = Z';

    % Distortion region (rotated ellipse)
    xc = Lx/2;
    zc = Lz/2;

    rx = 0.05 * Lx;   
    rz = 0.39 * Lz;

    theta = 38 * pi/180;
    
    R = [cos(theta)  sin(theta);
        -sin(theta)  cos(theta)];
    
    % Level of distortion is adjustable!
    strength = 0.3;
    
    for i = 2:nx
        for j = 2:nz
    
            x = X(i,j);
            z = Z(i,j);
    
            % shift to center
            v = [x-xc; z-zc];
    
            % rotate
            vr = R * v;
    
            xr = vr(1);
            zr = vr(2);
    
            % ellipse in rotated frame
            r = (xr/rx)^2 + (zr/rz)^2;
    
            if r < 1
                w = (1 - r);
    
                dx = strength * w * (x - xc);
                dz = strength * w * (z - zc);
    
                X(i,j) = x + dx;
                Z(i,j) = z + dz;
            end
        end
    end

    vertices = [X(:), Z(:)];

    % Initialize containers
    face_struct = struct('cells', {}, 'normal', {}, 'center', {}, 'area', {});
    face_map = containers.Map;
    face_counter = 0;

    cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'volume', {});
    cells = {};

    % Build cells
    cell_id = 0;

    for j = 1:nz
        for i = 1:nx
            cell_id = cell_id + 1;

            v1 = (j-1)*(nx+1) + i;
            v2 = v1 + 1;
            v3 = v2 + (nx+1);
            v4 = v1 + (nx+1);
            vids = [v1 v2 v3 v4];
            cells{cell_id} = vids;

            pts = vertices(vids, :);

            % centroid
            xc_cell = mean(pts,1)';
            cell_struct(cell_id).center = xc_cell;

            % area (shoelace)
            x = pts(:,1);
            z = pts(:,2);
            A = 0.5 * abs(sum(x .* circshift(z,-1)) - sum(z .* circshift(x,-1)));
            cell_struct(cell_id).volume = A;

            % Faces
            face_ids = zeros(1,4);
            face_signs = zeros(1,4);

            edge_list = {[v1 v2], [v2 v3], [v3 v4], [v4 v1]};

            for k = 1:4
                edge = edge_list{k};
                edge_sorted = sort(edge);
                key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));

                if isKey(face_map, key)
                    f = face_map(key);
                    face_struct(f).cells(end+1) = cell_id;
                else
                    face_counter = face_counter + 1;
                    f = face_counter;
                    face_map(key) = f;

                    p1 = vertices(edge_sorted(1),:);
                    p2 = vertices(edge_sorted(2),:);
                    t = p2 - p1;

                    normal = [-t(2); t(1)];
                    normal = normal / norm(normal);

                    center = 0.5*(p1+p2);
                    area = norm(t);

                    face_struct(f).cells = cell_id;
                    face_struct(f).normal = normal;
                    face_struct(f).center = center;
                    face_struct(f).area = area;
                end

                face_ids(k) = f;

                nf = face_struct(f).normal;
                xf = face_struct(f).center(:);
                s = sign(dot(nf, xf - xc_cell));
                if s == 0, s = 1; end
                face_signs(k) = s;
            end

            cell_struct(cell_id).faces = face_ids;
            cell_struct(cell_id).faces_orientation = face_signs;
        end
    end

    fprintf('Distorted structured mesh built:\n');
    fprintf('  %d vertices\n', size(vertices,1));
    fprintf('  %d cells\n', length(cell_struct));
    fprintf('  %d faces\n', length(face_struct));
end

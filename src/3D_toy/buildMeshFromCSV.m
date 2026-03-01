function [cell_struct, face_struct, vertices, cells] = buildMeshFromCSV(pointFile, polygonFile)
% Construct a 2D polygonal mesh from SPE11 style CSV geometry files

    % Read geometry data
    vertices = readmatrix(pointFile);

    % SNAP vertices near top boundary 
    % NOTE: Snap near-top vertices to y = Lz to avoid boundary detection errors
    % caused by small geometric perturbations in the input CSV mesh.
    Lz = 1200.02;     
    tol_snap = 1.0;
    
    idx = abs(vertices(:,2) - Lz) < tol_snap;
    fprintf('Snapping %d vertices to top boundary (z = %.5f)\n', sum(idx), Lz);
    vertices(idx,2) = Lz;
    raw_cells = readmatrix(polygonFile);

    % Remove NaNs and store vertex indices for each polygon
    nCells = size(raw_cells, 1);
    cells = cell(nCells, 1);
    for c = 1:nCells
        valid_ids = raw_cells(c, ~isnan(raw_cells(c, :)));
        cells{c} = valid_ids + 1;   % Convert 0-based to 1-based indexing
    end

    % Initialize structures
    cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'volume', {});
    face_struct = struct('cells', {}, 'normal', {}, 'center', {}, 'area', {});
    face_map = containers.Map;
    face_counter = 0;

    % Loop over cells
    for c = 1:nCells
        vids = cells{c};
        pts = vertices(vids, :);

        % Cell center
        xc = mean(pts, 1)';
        cell_struct(c).center = xc;

        % Polygon area & centroid (shoelace)
        x = pts(:,1);
        z = pts(:,2);

        x2 = circshift(x, -1);
        z2 = circshift(z, -1);

        cross = x .* z2 - x2 .* z;

        A = 0.5 * sum(cross);

        Cx = sum((x + x2) .* cross) / (6 * A);
        Cz = sum((z + z2) .* cross) / (6 * A);

        xc = [Cx; Cz];
        cell_struct(c).center = xc;
        cell_struct(c).volume = abs(A);

        % Build face list
        nverts = numel(vids);
        face_ids = zeros(1, nverts);
        face_signs = zeros(1, nverts);

        for k = 1:nverts
            v1 = vids(k);
            v2 = vids(mod(k, nverts) + 1);
            edge_sorted = sort([v1, v2]);
            key = sprintf('%d_%d', edge_sorted(1), edge_sorted(2));

            % If face already exists, append cell
            if isKey(face_map, key)
                f = face_map(key);
                face_struct(f).cells(end+1) = c;
            else
                % New face
                face_counter = face_counter + 1;
                f = face_counter;
                face_map(key) = f;

                p1 = vertices(edge_sorted(1), :);
                p2 = vertices(edge_sorted(2), :);
                t = p2 - p1;
                normal = [-t(2); t(1)];
                normal = normal / norm(normal);
                center = 0.5 * (p1 + p2);
                length_face = norm(t);

                face_struct(f).cells = c;
                face_struct(f).normal = normal;
                face_struct(f).center = center;
                face_struct(f).area = length_face;
            end

            face_ids(k) = f;

            % Orientation sign (dot product with face normal)
            nf = face_struct(f).normal;
            xf = face_struct(f).center(:);
            face_signs(k) = sign(dot(nf, xf - xc));
        end

        cell_struct(c).faces = face_ids;
        cell_struct(c).faces_orientation = face_signs;
    end

    fprintf('Mesh construction complete:\n');
    fprintf('  %d vertices\n', size(vertices,1));
    fprintf('  %d cells\n', nCells);
    fprintf('  %d faces\n', numel(face_struct));
end
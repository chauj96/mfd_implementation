function plotStreamlinesFromFlux(cell_struct, face_struct, flux)

    normals = cell2mat(arrayfun(@(f) f.normal(:), face_struct, 'UniformOutput', false))';
    areas = cell2mat(arrayfun(@(f) f.area, face_struct, 'UniformOutput', false))';
    xc = cell2mat(arrayfun(@(c) c.center(:), cell_struct, 'UniformOutput', false))';

    assert(size(normals,2) == 2, 'Only works for 2D grids.');
    assert(length(flux) == length([face_struct.normal]), 'Flux size must match number of faces.');


    % Convert normal flux to vectorial face fluxes
    v_faces = flux .* normals;  % [n_faces x 2]

    % Interpolate to cell centers by averaging neighboring face vectors
    nc = length(cell_struct);
    v_cells = zeros(nc, 2);
    face_counts = zeros(nc, 1);

    for c = 1:nc
        faces = cell_struct(c).faces;  % face indices
        for f = 1:length(faces)
            v_cells(c,:) = v_cells(c,:) + v_faces(faces(f), :);
        end
        face_counts(c) = length(faces);
    end

    v_cells = v_cells ./ face_counts;  % average

    % Interpolate to regular grid for streamslice
    x = xc(:,1); 
    y = xc(:,2);
    vx = v_cells(:,1); 
    vy = v_cells(:,2);

    [xq, yq] = meshgrid(linspace(min(x), max(x), nc), ...
                        linspace(min(y), max(y), nc));

    % Interpolate vector field
    vq_x = griddata(x, y, vx, xq, yq, 'linear');
    vq_y = griddata(x, y, vy, xq, yq, 'linear');

    % Plot streamlines
    figure;
    streamslice(xq, yq, vq_x, vq_y);
    hold on;
    %quiver(x, y, vx, vy, 'r');  % optional: vector field at cells
    axis equal tight;
    title('2D Streamlines from Normal Facet Fluxes');
    legend('Streamlines');
end

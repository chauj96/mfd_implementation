function [q_faces, p_cells] = projectExactSolution(cell_struct, face_struct)
    % Evaluate exact pressures at cell centers and fluxes across faces
    %
    % p_exact(x, y) = 1 - x
    % flux_exact = [1; 0] (constant vector field)
    
    n_cells = numel(cell_struct);
    n_faces = numel(face_struct);

    p_cells = zeros(n_cells, 1);
    q_faces = zeros(n_faces, 1);

    % Evaluate pressure at cell centers
    for c = 1:n_cells
        xc = cell_struct(c).center(1);  % x-coordinate
        p_cells(c) = 1 - xc;
    end

    % Evaluate flux across faces
    for f = 1:n_faces
        n = face_struct(f).normal(:);      % face normal (column vector)
        a = face_struct(f).area;
        flux_exact = [1; 0];               % constant flux
        q_faces(f) = dot(flux_exact, n * a);
    end
end

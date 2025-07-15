function R = compute_residual(Sw_new, Sw_old, m, phi, dt, dx)
    n_cells = length(Sw_new);
    n_faces = n_cells + 1;

    % Upwinding for f_w at faces
    fw_faces = zeros(n_faces, 1);

    % Left boundary face (i = 1)
    if m(1) >= 0
        fw_faces(1) = 1.0; % inject pure water
    else
        fw_faces(1) = Sw_new(1); % inflow from cell
    end

    % Interior faces (i = 2 to n_faces - 1)
    for i = 2:n_faces-1
        if m(i) >= 0
            fw_faces(i) = Sw_new(i-1);
        else
            fw_faces(i) = Sw_new(i);
        end
    end

    % Right boundary face (i = n_faces)
    if m(1) >= 0
        fw_faces(end) = Sw_new(end); % outflow from last cell
    else
        fw_faces(end) = 0.0;
    end

    % Compute flux term
    flux = fw_faces .* m;

    % Residual vector (n_cells * 1)
    R = phi * (Sw_new - Sw_old) * dx / dt + (flux(2:end) - flux(1:end-1));
end
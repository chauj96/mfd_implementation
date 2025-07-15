function [cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct)

    % constants
    K_vals = 1.0; % permeability
    phi_vals = 0; % porosity
    rho_vals = 1000; % fluid density [kg/m^3]
    g_val = 0.0; % gravitational acceleration [m/s^2]
    gravity_dir = [0; -1];

    % fetch face centers
    f_centers = vertcat(face_struct.center);

    bottom_idx = find(f_centers(:,2) == 0.0);
    top_idx = find(f_centers(:,2) == 2.0);
    BC_Dirichlet_map = containers.Map([top_idx, bottom_idx], [repmat(4, 1, length(top_idx)),repmat(2, 1, length(bottom_idx))]);

    east_idx = find(f_centers(:,1) == 0.0);
    west_idx = find(f_centers(:,1) == 1.0);
    BC_Neumann_map = containers.Map([east_idx, west_idx], [repmat(0.0, 1, length(east_idx)),repmat(0.0, 1, length(west_idx))]);
    
    % assign to each cell
    for c = 1:length(cell_struct)
        cell_struct(c).K = K_vals;
        cell_struct(c).phi = phi_vals;
        cell_struct(c).rho = rho_vals;
    end
    
    % assign gravity and density to each face
    for f = 1:length(face_struct)

        face_struct(f).gravity = g_val * gravity_dir;
        face_struct(f).rho = rho_vals;

        if isKey(BC_Dirichlet_map, f)
            face_struct(f).BC_pressure = BC_Dirichlet_map(f);
        end

        if isKey(BC_Neumann_map, f)
            face_struct(f).BC_flux = BC_Neumann_map(f);
        end
        

    end

end
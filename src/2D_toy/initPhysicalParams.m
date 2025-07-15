function [cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct)

    % constants
    K_vals = 1.0; % permeability
    phi_vals = 0; % porosity
    rho_vals = 1000; % fluid density [kg/m^3]
    g_val = 0.0; % gravitational acceleration [m/s^2]
    gravity_dir = [0; -1];
    % Later change to "map"
    BC_data = [4, 0.5, 2.0; 2, 0.5, 0.0];
    BC_coord = BC_data(:,2:end);

    BC_data_neumann = [0, 0.0, 0.5; 0, 1.0, 0.5; 0, 0.0, 1.5; 0, 1.0, 1.5];
    BC_coord_neumann = BC_data_neumann(:,2:end);
    
    % assign to each cell
    for c = 1:length(cell_struct)
        cell_struct(c).K = K_vals;
        cell_struct(c).phi = phi_vals;
        cell_struct(c).rho = rho_vals;
    end
    
    % assign gravity and density to each face
    for f = 1:length(face_struct)
        xc = face_struct(f).center;
        x_diff_D = [xc;xc] - BC_coord;
        x_diff_N = [xc;xc;xc;xc] - BC_coord_neumann;
        
        is_Member_D = vecnorm(x_diff_D') <= 1e-15;
        is_Member_N = vecnorm(x_diff_N') <= 1e-15;

        face_struct(f).gravity = g_val * gravity_dir;
        face_struct(f).rho = rho_vals;

        if any(is_Member_D, 'all')
            idx = find(is_Member_D == 1);
            face_struct(f).BC_pressure = BC_data(idx,1);
        else
            face_struct(f).BC_pressure = 0.0; 
        end

        if any(is_Member_N, 'all')
            idx = find(is_Member_N == 1);
            face_struct(f).BC_flux = BC_data_neumann(idx,1);
        end
        

    end

end
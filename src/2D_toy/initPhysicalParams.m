function [cell_struct, face_struct] = initPhysicalParams(cell_struct, face_struct, Lx, Lz, case_type)

    % constants
    K_tensor = [1e-3, 0.0;
                0.0, 5.0]; % permeability tensor
    phi_vals = 0; % porosity
    rho_vals = 1000; % fluid density [kg/m^3]
    g_val = 0.0; % gravitational acceleration [m/s^2]
    gravity_dir = [0; -1];
    tol = 1e-8;

    % fetch face centers
    if strcmp(case_type, 'structured')
        f_centers = vertcat(face_struct.center);
    elseif strcmp(case_type, 'unstructured')
        n_faces = length(face_struct);
        f_centers = zeros(n_faces,2) ; % preallocate
    
        for f = 1:n_faces
            f_centers(f, :) = face_struct(f).center(:)';
        end
    end

    bottom_idx = find(abs(f_centers(:,2) - 0.0) < tol);
    top_idx = find(abs(f_centers(:,2) - Lz) < tol);
    BC_Dirichlet_map = containers.Map([top_idx; bottom_idx], [repmat(4, 1, length(top_idx)),repmat(2, 1, length(bottom_idx))]);
 %   BC_Dirichlet_map = containers.Map([top_idx], [repmat(4, 1, length(top_idx))]);

    east_idx = find(abs(f_centers(:,1) - 0.0) < tol);
    west_idx = find(abs(f_centers(:,1) - Lx) < tol);
    BC_Neumann_map = containers.Map([east_idx; west_idx], [repmat(0.0, 1, length(east_idx)),repmat(0.0, 1, length(west_idx))]);
%    BC_Neumann_map = containers.Map([east_idx; west_idx; bottom_idx], [repmat(0.0, 1, length(east_idx)),repmat(0.0, 1, length(west_idx)),repmat(0.0, 1, length(bottom_idx))]);
    
    % assign to each cell
    for c = 1:length(cell_struct)
        cell_struct(c).K = K_tensor;

%         x = cell_struct(c).center(1);
%         if x < 0.1
%             cell_struct(c).K = [0.01, 0.01; 0.01, 0.01];
%         else
%             cell_struct(c).K = [10.0, 0.0; 0.0, 10.0];
%         end

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
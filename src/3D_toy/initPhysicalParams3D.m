function [cell_struct, face_struct, phys] = initPhysicalParams3D(cell_struct, face_struct, Lx, Ly, Lz, bc_option)

    % Physical constants
    K_tensor = eye(3);    
    phi_vals = 0.3;
    rho_vals = 1000;
    g_val = 0.0;
    gravity_dir = [0; 0; -1];
    tol = 1e-6;

    % Reference pressure gradient
    % (x-direction flow)
    % grad_pref = [-1/Lx; 0; 0];  
    % grad_pref = [-1/Lx; -1/Ly; -1/Lz];
    % m_ref_vec = -K_tensor * grad_pref;

    % Face centers
    n_faces = length(face_struct);
    f_centers = zeros(n_faces,3);
    for f = 1:n_faces
        f_centers(f,:) = face_struct(f).center(:)';
    end

    % Identify boundary faces
    west_idx = find(abs(f_centers(:,1) - 0.0) < tol);
    east_idx = find(abs(f_centers(:,1) - Lx ) < tol);

    south_idx = find(abs(f_centers(:,2) - 0.0) < tol);
    north_idx = find(abs(f_centers(:,2) - Ly ) < tol);

    bottom_idx = find(abs(f_centers(:,3) - 0.0) < tol);
    top_idx = find(abs(f_centers(:,3) - Lz ) < tol);

    %%
    if strcmp(bc_option, 'linear')
        % Reference pressure gradient
        % (x-direction flow)
        grad_pref = [-1/Lx; 0; 0];  
        m_ref_vec = -K_tensor * grad_pref;

        % Boundary conditions
        % Dirichlet: x = 0 -> p = 1, x = Lx -> p = 0
        BC_Dirichlet_map = containers.Map( ...
            [west_idx; east_idx], ...
            [ones(length(west_idx),1); zeros(length(east_idx),1)] );
    
        % Neumann: no-flow everywhere else
        BC_Neumann_map = containers.Map( ...
            [south_idx; north_idx; bottom_idx; top_idx], ...
            zeros(length([south_idx; north_idx; bottom_idx; top_idx]),1) );
        
    elseif strcmp(bc_option, 'corner2corner')
        % Reference pressure gradient
        grad_pref = [-1/Lx; -1/Ly; -1/Lz];
        m_ref_vec = -K_tensor * grad_pref;

        % Boundary faces (all exterior faces)
        boundary_faces = unique([ ...
            west_idx;
            east_idx;
            south_idx;
            north_idx;
            bottom_idx;
            top_idx ]);
    
        n_cells = length(cell_struct);
        cell_centers = zeros(n_cells,3);
    
        for c = 1:n_cells
            cell_centers(c,:) = cell_struct(c).center(:)';
        end
    
        [~, inlet_cell] = min(vecnorm(cell_centers - [0, 0, Lz], 2, 2));
        [~, outlet_cell] = min(vecnorm(cell_centers - [Lx, Ly, 0], 2, 2));
    
        % Dirichlet faces = exposed boundary faces of those cells
        inlet_faces_all  = cell_struct(inlet_cell).faces(:);
        outlet_faces_all = cell_struct(outlet_cell).faces(:);
    
        inlet_faces  = inlet_faces_all(ismember(inlet_faces_all, boundary_faces));
        outlet_faces = outlet_faces_all(ismember(outlet_faces_all, boundary_faces));
    
        dirichlet_faces = unique([inlet_faces; outlet_faces]);
    
        % Pressure values from analytical solution p=ax+by+cz+d
        BC_Dirichlet_vals = zeros(length(dirichlet_faces),1);
    
        for k = 1:length(dirichlet_faces)
            f = dirichlet_faces(k);
    
            xf = face_struct(f).center(1);
            yf = face_struct(f).center(2);
            zf = face_struct(f).center(3);
    
            BC_Dirichlet_vals(k) = ...
                grad_pref(1)*xf + ...
                grad_pref(2)*yf + ...
                grad_pref(3)*zf + 1;
        end
    
        BC_Dirichlet_map = containers.Map( ...
            dirichlet_faces, ...
            BC_Dirichlet_vals );
    
        % Remaining boundary faces = no-flow Neumann
        neumann_faces = setdiff(boundary_faces, dirichlet_faces);
    
        BC_Neumann_map = containers.Map( ...
            neumann_faces, ...
            zeros(length(neumann_faces),1) );
    end

    %%

    % Assign cell properties
    for c = 1:length(cell_struct)
        cell_struct(c).K = K_tensor;
        cell_struct(c).phi = phi_vals;
        cell_struct(c).rho = rho_vals;
    end

    for f = 1:n_faces
        face_struct(f).BC_flux = [];
    end

    % Assign face properties + BC
    for f = 1:n_faces
        face_struct(f).gravity = g_val * gravity_dir;
        face_struct(f).rho = rho_vals;

        if isKey(BC_Dirichlet_map, f)
            face_struct(f).BC_pressure = BC_Dirichlet_map(f);
        end

        if isKey(BC_Neumann_map, f)
            face_struct(f).BC_flux = BC_Neumann_map(f);
        end
    end

    % Analytical flux on faces
    m_ref_faces = zeros(n_faces,1);
    for f = 1:n_faces
        n_f = face_struct(f).normal(:);
        n_f = n_f / norm(n_f);
        Af = face_struct(f).area;

        m_ref_faces(f) = Af * dot(m_ref_vec, n_f);
    end

    % Store phys struct (actually m_ref_faces is the same as m_proj)
    phys.K_tensor = K_tensor;
    phys.grad_pref = grad_pref;
    phys.m_ref_vec = m_ref_vec;
    phys.m_ref_faces = m_ref_faces;
end
function [cell_struct, face_struct, phys] = initPhysicalParams3D(cell_struct, face_struct, Lx, Ly, Lz)

    % Physical constants
    K_tensor = eye(3);    
    phi_vals = 0.0;
    rho_vals = 1000;
    g_val = 0.0;
    gravity_dir = [0; 0; -1];
    tol = 1e-6;

    % Reference pressure gradient
    % (x-direction flow)
    grad_pref = [-1/Lx; 0; 0];         
    m_ref_vec = -K_tensor * grad_pref;

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

    % Boundary conditions
    % Dirichlet: x = 0 -> p = 1, x = Lx -> p = 0
    BC_Dirichlet_map = containers.Map( ...
        [west_idx; east_idx], ...
        [ones(length(west_idx),1); zeros(length(east_idx),1)] );

    % Neumann: no-flow everywhere else
    BC_Neumann_map = containers.Map( ...
        [south_idx; north_idx; bottom_idx; top_idx], ...
        zeros(length([south_idx; north_idx; bottom_idx; top_idx]),1) );

    % Assign cell properties
    for c = 1:length(cell_struct)
        cell_struct(c).K = K_tensor;
        cell_struct(c).phi = phi_vals;
        cell_struct(c).rho = rho_vals;
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
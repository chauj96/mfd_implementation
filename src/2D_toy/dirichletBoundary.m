function rhs_Dirichlet = dirichletBoundary(cell_struct, face_struct)

    n_cells = length(cell_struct);
    n_faces = length(face_struct);
    
    rhs_Dirichlet = zeros(n_faces,1);
    
    for k = 1:n_cells
        face_ids = cell_struct(k).faces;
        signs = cell_struct(k).faces_orientation;
    
        for j = 1:length(face_ids)
            
            f = face_ids(j);
            if isempty(face_struct(f).BC_pressure)
                continue;
            end
    
            sigma = signs(j);
            area_f = face_struct(f).area;
            p_D = face_struct(f).BC_pressure;
            rhs_Dirichlet(f) = sigma * p_D * area_f;
        end
    end

end
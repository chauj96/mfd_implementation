function [A, rhs_BC] = neumannBoundary(A, cell_struct, face_struct)

    n_faces = length(face_struct);
    rhs_BC = zeros(n_faces,1);
    
    f_ids = find(arrayfun(@(x) ~isempty(x.BC_flux), face_struct));
    
    for idx = 1:length(f_ids)
        f = f_ids(idx);
        A(f, :) = 0;
        A(:, f) = 0;
        A(f, f) = 1;
    
        rhs_BC(f) = face_struct(f).BC_flux;
    end
end
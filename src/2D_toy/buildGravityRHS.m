function f_g = buildGravityRHS(face_struct, g)
    n_faces = length(face_struct);
    f_g = zeros(n_faces, 1);
    
    for f = 1:n_faces
        rho = face_struct(f).rho;
        n = face_struct(f).normal(:);
        area_f = face_struct(f).area;
        f_g(f) = -rho * dot([0; -g], n) * area_f;
    end
end
function f_g = buildGravityRHS(face_struct, g)
n_faces = length(face_struct);
f_g = zeros(n_faces, 1);

for f = 1:n_faces
    rho = face_struct(f).rho;
    n = face_struct(f).normal(:);
    f_g(f) = face_struct(f).invT * rho * dot([0; -g], n);
end

end
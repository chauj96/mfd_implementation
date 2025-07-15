function BC_vec = dirichletBoundary(cell_struct, face_struct)

n_cells = length(cell_struct);
n_faces = length(face_struct);

BC_vec = zeros(n_faces,1);

for k = 1:n_cells
    face_ids = cell_struct(k).faces;
    signs = cell_struct(k).face_dirs;

    for j = 1:length(face_ids)
        f = face_ids(j);
        sigma = signs(j);
        area_f = face_struct(f).area;
        pressure = face_struct(f).BC_pressure;

        BC_vec(f) = sigma * area_f * pressure;
    end
end

end
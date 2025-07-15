function B = buildBmatrix(cell_struct, face_struct)

n_cells = length(cell_struct);
n_faces = length(face_struct);

rows = [];
cols = [];
vals = [];

for k = 1:n_cells
    face_ids = cell_struct(k).faces;
    signs = cell_struct(k).face_dirs;

    for j = 1:length(face_ids)
        f = face_ids(j);
        sigma = signs(j);
        area_f = face_struct(f).area;

        rows(end+1) = k;
        cols(end+1) = f;
        vals(end+1) = sigma * area_f;
    end
end

B = sparse(rows, cols, vals, n_cells, n_faces);
end
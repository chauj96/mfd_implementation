function M = buildMmatrixMFEM(cell_struct,face_struct)
% Build the M matrix (inner product for velocity DOFs)

rows = [];
cols = [];
vals = [];

n_cells = length(cell_struct);
n_faces = length(face_struct);

for c = 1:n_cells
    face_ids = cell_struct(c).faces;
    for j = 1:length(face_ids)

        f = face_ids(j);
        xC = cell_struct(c).center(:);
        xF = face_struct(f).center(:);
        nf = face_struct(f).normal(:); 
        Af = face_struct(f).area;  
        K = cell_struct(c).K;
        d = abs(dot(xF - xC, nf));
        invT = d / (K * Af);

        rows(end+1) = f;
        cols(end+1) = f;
        vals(end+1) = invT;
    end
end

M = sparse(rows, cols, vals, n_faces, n_faces);

end
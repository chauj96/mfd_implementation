function B = buildBmatrix_1d(cell_struct)
% Build the B matrix for 1D fracture (discrete divergence operator)

    n_cells = length(cell_struct);
    n_faces = n_cells + 1;

    % Precompute total number of non-zeros
    total_nnz = sum(cellfun(@(x) length(x), {cell_struct.faces}));

    rows = zeros(total_nnz, 1);
    cols = zeros(total_nnz, 1);
    vals = zeros(total_nnz, 1);

    idx = 0;
    for k = 1:n_cells
        face_ids = cell_struct(k).faces;
        global_face_ids = [(k-1)+1, (k-1)+2];
        signs = cell_struct(k).faces_orientation;

        for j = 1:length(face_ids)
            idx = idx + 1;
            rows(idx) = k;
            cols(idx) = global_face_ids(j);
            vals(idx) = signs(j);
        end
    end

    B = sparse(rows, cols, vals, n_cells, n_faces);
end
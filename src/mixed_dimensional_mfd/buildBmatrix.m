function B = buildBmatrix(cell_struct, face_struct)
    % Build the B matrix (discrete divergence operator)

    n_cells = length(cell_struct);
    n_faces = length(face_struct);

    % Precompute total number of non-zeros
    total_nnz = sum(cellfun(@(x) length(x), {cell_struct.faces}));

    rows = zeros(total_nnz, 1);
    cols = zeros(total_nnz, 1);
    vals = zeros(total_nnz, 1);

    idx = 0;
    for k = 1:n_cells
        face_ids = cell_struct(k).faces;
        signs = cell_struct(k).faces_orientation;

        for j = 1:length(face_ids)
            idx = idx + 1;
            rows(idx) = k;
            cols(idx) = face_ids(j);
            vals(idx) = signs(j);
        end
    end

    B = sparse(rows, cols, vals, n_cells, n_faces);
end

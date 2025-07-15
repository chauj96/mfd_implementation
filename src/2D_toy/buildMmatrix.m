function M = buildMmatrix(face_struct)
% Build the M matrix (inner product for velocity DOFs)
% M: diagonal matrix with T_f values

n_faces = length(face_struct);
invT_vec = zeros(n_faces,1);

for i = 1:n_faces
    invT_vec(i) = face_struct(i).invT;
end

M = spdiags(invT_vec, 0, n_faces, n_faces);
end
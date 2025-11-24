function M = applyCoupling(M, cell_struct_2d, cell_struct_1d, d)

    idx = size(M,1) - (length(cell_struct_1d) - 1);
    for i = 1:length(cell_struct_1d)
       left_cell = cell_struct_1d(i).parent_2d(1);
       right_cell = cell_struct_1d(i).parent_2d(2);

       face1 = cell_struct_2d(left_cell).faces(2); % select right face
       face2 = cell_struct_2d(right_cell).faces(4); % select left face

       % Add T and T^t matrices 
       M(idx, face1) = -1;  % left face of current fracture
       M(idx, face2) = 1;   % right face of current fracture
       M(face1, idx) = -1;
       M(face2, idx) = 1;

       % Add Robin coupling to M_2d matrix: alpha^{-1} = d / (2 * K_ortho)
       M(face1, face1) = M(face1, face1) + d / (2 * cell_struct_1d(i).K(1,1));
       M(face2, face2) = M(face2, face2) + d / (2 * cell_struct_1d(i).K(1,1));

       idx = idx + 1;
    end

end
function cell_struct = createBmatrix(cell_struct)

    for k = 1:length(cell_struct)
        cell_struct(k).B = cell_struct(k).faces_orientation(:);  
    end

end
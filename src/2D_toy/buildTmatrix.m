function T = buildTmatrix(cell_struct)
% T: sparse diagonal matrix (n_cells * n_cells)

n_cells = length(cell_struct);
t_vec = zeros(n_cells,1);

for i = 1:n_cells
    phi = cell_struct(i).phi; % porosity
    V = cell_struct(i).volume; % cell volume
    t_vec(i) = phi * V;
end

T = spdiags(t_vec, 0, n_cells, n_cells);

end
function plotPressurePolygonal(vertices, cells, p_sol)
    n_cells = length(cells);

    figure;
    hold on;

    for c = 1:n_cells
        vert_ids = cells{c};
        poly = vertices(vert_ids, :);

        fill(poly(:,1), poly(:,2), p_sol(c), 'EdgeColor', 'k');
    end
    
    axis equal tight;
    colorbar;
    title('Pressure Field');
    xlabel('x');
    ylabel('z');
end

% No longer use this plot function (move to plotPressurePolygonal.m)

function plotPressure(cell_struct, p_sol)
    % Extract grid info
    n_cells = length(cell_struct);
    xc = zeros(n_cells, 1);
    zc = zeros(n_cells, 1);

    for i = 1:n_cells
        xc(i) = cell_struct(i).center(1);
        zc(i) = cell_struct(i).center(2);
    end

    % Scatter plot (colored by pressure)
    scatter(xc, zc, 100, p_sol, 'filled');
%     plot(zc, p_sol, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    colorbar;
    xlabel('z'); ylabel('pressue');
    title('Pressure field');
    axis equal tight;
end
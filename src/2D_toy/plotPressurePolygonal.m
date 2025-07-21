function plotPressurePolygonal(vertices, cells, p_sol)
    n_cells = length(cells);

    figure;
    hold on;

    % Plot filled polygonal cells
    cell_centers = zeros(n_cells, 2);
    for c = 1:n_cells
        vert_ids = cells{c};
        poly = vertices(vert_ids, :);
        fill(poly(:,1), poly(:,2), p_sol(c), 'EdgeColor', [0.7, 0.7, 0.7]);
        cell_centers(c,:) = mean(poly, 1); % Store cell centroid
    end
    
    % Set axis
    axis equal tight;

    % Build interpolation grid
    x_min = min(vertices(:,1));
    x_max = max(vertices(:,1));
    y_min = min(vertices(:,2));
    y_max = max(vertices(:,2));
    
    [X, Y] = meshgrid( ...
        linspace(x_min, x_max, 300), ...
        linspace(y_min, y_max, 300));
    
    % Interpolate the pressure field to the grid
    F = scatteredInterpolant(cell_centers(:,1), cell_centers(:,2), p_sol(:), 'natural', 'none');
    P_interp = F(X, Y);

    % Plot colored and thick contour lines
    hold on;
    [C, h] = contour(X, Y, P_interp, 13, 'LineColor', [1.0, 0.0, 0.0]);
    h.LineWidth = 3;                     % Thicker lines
    
    clabel(C, h, 'FontSize', 12, 'Color', 'k');

    % Use the same colormap for contour and fill
    colormap(parula);
    colorbar;

    title('Pressure Field with Contours (Polygonal Mesh)');
    xlabel('x');
    ylabel('z');
end

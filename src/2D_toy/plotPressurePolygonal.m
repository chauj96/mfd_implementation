function plotPressurePolygonal(vertices, cells, p_sol, title_suffix, scheme_labels)
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

    % Add scheme labels ("T" or "M") on top of each cell
    if nargin >= 5
        for c = 1:n_cells
            if scheme_labels{c} == 'T'
                color = [0.5 0.5 0.5]; % blue for TPFA
            else
                color = [0.5 0.5 0.5];  % red for MFD
            end
            text(cell_centers(c,1), cell_centers(c,2), scheme_labels{c}, ...
                'HorizontalAlignment','center', 'FontSize',12, ...
                'FontWeight','bold', 'Color', color);
        end
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
    F = scatteredInterpolant(cell_centers(:,1), cell_centers(:,2), p_sol(:), 'linear', 'none');
    P_interp = F(X, Y);

    % Plot colored and thick contour lines
    hold on;
    [C, h] = contour(X, Y, P_interp, 9, 'LineColor', [1.0, 0.0, 0.0]);
    h.LineWidth = 3;
    clabel(C, h, 'FontSize', 12, 'Color', 'k');

    % Apply red-blue colormap
    colormap(gca, redblue(9));
    caxis([0 1]);
    colorbar;

    % Title and labels
    title(['Pressure Field: ', title_suffix]);
    xlabel('x'); ylabel('z');

    %% Helper: redblue colormap function
    function cmap = redblue(m)
        if nargin < 1
            m = size(get(gcf, 'colormap'), 1);
        end
        bottom = [0 0 1]; middle = [1 1 1]; top = [1 0 0];
        n_half = floor(m/2);
        if mod(m,2)==0
            r = [linspace(bottom(1), middle(1), n_half), linspace(middle(1), top(1), n_half)];
            g = [linspace(bottom(2), middle(2), n_half), linspace(middle(2), top(2), n_half)];
            b = [linspace(bottom(3), middle(3), n_half), linspace(middle(3), top(3), n_half)];
        else
            r = [linspace(bottom(1), middle(1), n_half), middle(1), linspace(middle(1), top(1), n_half)];
            g = [linspace(bottom(2), middle(2), n_half), middle(2), linspace(middle(2), top(2), n_half)];
            b = [linspace(bottom(3), middle(3), n_half), middle(3), linspace(middle(3), top(3), n_half)];
        end
        cmap = [r', g', b'];
    end
end

% function plotPressurePolygonal(vertices, cells, p_sol, title_suffix)
%     n_cells = length(cells);
% 
%     figure;
%     hold on;
% 
%     % Plot filled polygonal cells
%     cell_centers = zeros(n_cells, 2);
%     for c = 1:n_cells
%         vert_ids = cells{c};
%         poly = vertices(vert_ids, :);
%         fill(poly(:,1), poly(:,2), p_sol(c), 'EdgeColor', [0.7, 0.7, 0.7]);
%         cell_centers(c,:) = mean(poly, 1); % Store cell centroid
%     end
% 
%     % Set axis
%     axis equal tight;
% 
%     % Build interpolation grid
%     x_min = min(vertices(:,1));
%     x_max = max(vertices(:,1));
%     y_min = min(vertices(:,2));
%     y_max = max(vertices(:,2));
% 
%     [X, Y] = meshgrid( ...
%         linspace(x_min, x_max, 300), ...
%         linspace(y_min, y_max, 300));
% 
%     % Interpolate the pressure field to the grid
%     F = scatteredInterpolant(cell_centers(:,1), cell_centers(:,2), p_sol(:), 'linear', 'none');
%     P_interp = F(X, Y);
% 
%     % Plot colored and thick contour lines
%     hold on;
%     [C, h] = contour(X, Y, P_interp, 9, 'LineColor', [1.0, 0.0, 0.0]);
%     h.LineWidth = 3;                     % Thicker lines
% 
%     clabel(C, h, 'FontSize', 12, 'Color', 'k');
% 
%     % Use the same colormap for contour and fill
%     % Apply custom red-blue colormap
%     colormap(gca, redblue(9)); % Call the custom function defined below
%     caxis([0 1]);
%     colorbar;
% 
%     %% === Helper: redblue colormap function ===
%     function cmap = redblue(m)
%         if nargin < 1
%             m = size(get(gcf, 'colormap'), 1);
%         end
%         bottom = [0 0 1];       % blue
%         middle = [1 1 1];       % white
%         top    = [1 0 0];       % red
% 
%         % Interpolation for smooth transition
%         n_half = floor(m/2);
%         if mod(m,2)==0
%             % Even number of levels
%             r = [linspace(bottom(1), middle(1), n_half), linspace(middle(1), top(1), n_half)];
%             g = [linspace(bottom(2), middle(2), n_half), linspace(middle(2), top(2), n_half)];
%             b = [linspace(bottom(3), middle(3), n_half), linspace(middle(3), top(3), n_half)];
%         else
%             % Odd number of levels
%             r = [linspace(bottom(1), middle(1), n_half), middle(1), linspace(middle(1), top(1), n_half)];
%             g = [linspace(bottom(2), middle(2), n_half), middle(2), linspace(middle(2), top(2), n_half)];
%             b = [linspace(bottom(3), middle(3), n_half), middle(3), linspace(middle(3), top(3), n_half)];
%         end
%         cmap = [r', g', b'];
%     end
% 
%     title(['Pressure Field: ', title_suffix]);
%     xlabel('x');
%     ylabel('z');
% end

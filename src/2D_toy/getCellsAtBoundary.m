function [cell_ids, z_vals] = getCellsAtBoundary(cell_struct, location)
% Return cell indices at top or bottom boundary

z_vals = arrayfun(@(c) c.center(2), cell_struct);
tol = 1e-8;

switch lower(location)
    case 'top'
        z_top = max(z_vals);
        cell_ids = find(abs(z_vals - z_top) < tol);
    case 'bottom'
        z_bot = min(z_vals);
        cell_ids = find(abs(z_vals - z_bot) < tol);
    otherwise
        error('Unknown boundary location');
end


end
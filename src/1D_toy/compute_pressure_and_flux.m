function [p_exact, p, m] = compute_pressure_and_flux(n_cells, dx, plot_flag)
    if nargin < 3
        plot_flag = false;
    end

    % Grid setup
    n_faces = n_cells + 1;
    p_left = 1.0;
    p_right = 0.0;

    % Operators
    G = zeros(n_faces, n_cells);
    for i = 1:n_faces
        if i > 1
            G(i, i-1) = -1;
        end
        if i <= n_cells
            G(i, i) = G(i, i) +1;
        end
    end
    G = G / dx;
    

    M = eye(n_faces);
    M_inv = inv(M);
    D = -G';

    % Boundary vector
    b_bc = zeros(n_faces, 1);
    b_bc(1) = -p_left / dx;
    b_bc(end) = p_right / dx;

    % Solve for pressure
    A = D * M_inv * G;
    b = -D * M_inv * b_bc;


    p = A \ b;

    % Flux at faces
    m = -M_inv * (G * p + b_bc);
    

    % Optional plot
    if plot_flag
        x_centers = linspace(dx/2, 1 - dx/2, n_cells)';
        p_exact = p_left - (p_left - p_right) * x_centers;
        figure;
        plot(x_centers, p, 'o-', 'DisplayName', 'MFD pressure');
        hold on;
        plot(x_centers, p_exact, 'k--', 'DisplayName', 'Exact pressure');
        xlabel('x'); ylabel('Pressure');
        title('Pressure profile: MFD vs Exact');
        legend('Location','best');
        grid on;
    end
end

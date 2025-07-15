% This is where we can run the whole simulation (1D toy problem)
clear; clc; close all;

% Parameters
n_cells = 2000;
L = 1.0;
dx = L / n_cells;
dt = 0.01;
nt = 200;
phi = 1.0;
tol = 1e-7;
max_iter = 20;

% Call pressure and flux values
[p_exact, p, m] = compute_pressure_and_flux(n_cells, dx, true);

% Initial saturation
Sw = zeros(n_cells, 1);
Sw(1) = 1.0;
Sw_hist = Sw';
mass_hist = phi * dx * sum(Sw);

% Time stepping
for t = 1:nt
    Sw_old = Sw;
    Sw_new = Sw_old;  % Initial guess

    % Newton iteration
    for iter = 1:max_iter
        R = compute_residual(Sw_new, Sw_old, m, phi, dt, dx);
        J = compute_jacobian(Sw_new, m, phi, dt, dx);

        delta = -J \ R;
        Sw_new = Sw_new + delta;

        % Check convergence
        if norm(delta, inf) < tol
            break;
        end
    end

    if iter == max_iter
        warning('Newton did not converge at timestep %d', t);
    end

    Sw = Sw_new;
    Sw_hist = [Sw_hist; Sw'];
    mass_hist = [mass_hist; phi * dx * sum(Sw)];
end

% Saturation plot
x = linspace(dx/2, L - dx/2, n_cells);
figure;
hold on;
plot(x, Sw_hist(1,:), 'k-', 'DisplayName', 't=0.00', 'LineWidth', 1);
for t = 1:5:nt+1
    plot(x, Sw_hist(t,:), '-', 'DisplayName', sprintf('t=%.2f', (t-1)*dt));
end
xlabel('x'); ylabel('S_w');
title('Change of Water Saturation over time');
legend('Location','bestoutside');
legend show;
grid on;

% Total water flux plot
time = linspace(0, nt*dt, nt+1);
figure;
plot(time, mass_hist, 'bo-', 'LineWidth', 1.0);
xlabel('Time');
ylabel('Total Water Mass');
title('Total Water Mass over time');
grid on;





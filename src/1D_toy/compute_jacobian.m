function J = compute_jacobian(Sw, m, phi, dt, dx)
    n = length(Sw);
    J = zeros(n,n);

    fw_prime = @(S) 1.0; % now we have fw = S => dfw/dS = 1

    for i = 1:n
        % Diagonal element: Accumulation term
        J(i,i) = phi * dx / dt;

        % Upwinding for flux derivative
        % Left face
        if i > 1
            if m(i) >= 0
                % Upwind: S_{i-1}
                J(i, i-1) = J(i, i-1) - m(i) * fw_prime(Sw(i-1));
            else
                % Upwind: S_{i}
                J(i, i) = J(i, i) - m(i) * fw_prime(Sw(i));
            end
        else
        end

        % Right face
        if i < n
            if m(i+1) >= 0
                % Upwind: S_{i}
                J(i, i) = J(i, i) + m(i+1) * fw_prime(Sw(i));
            else
                % Upwind: S_{i+1}
                J(i, i+1) = J(i, i+1) + m(i+1) * fw_prime(Sw(i+1));
            end
        else
            % Right boundary (outflow): Upwind S_{n}
            J(i, i) = J(i, i) + m(i+1) * fw_prime(Sw(i));
        end
    end
end
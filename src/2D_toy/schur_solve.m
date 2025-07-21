function sol = schur_solve(J, rhs, primary_idx, secondary_idx)
% with J = [A, B; C, D] using Schur complement elimination of secondary variables.
%
% INPUTS:
%   A             - full block matrix
%   rhs           - right-hand side vector
%   primary_idx   - indices for primary variables (e.g., cells)
%   secondary_idx - indices for secondary variables (e.g., faces)
%
% OUTPUT:
%   sol - full solution vector such that A * sol = rhs

    % Split matrix blocks
    A = J(secondary_idx, secondary_idx);
    B = J(secondary_idx, primary_idx);
    C = J(primary_idx, secondary_idx);
    D = J(primary_idx, primary_idx);

    % Split RHS
    f_s = rhs(secondary_idx);
    f_p = rhs(primary_idx);

    % Solve for secondary component
    Hinv_f_s = A \ f_s;
    Hinv_B   = A \ B;

    % Schur complement and reduced RHS
    S = D - C * Hinv_B;
    rhs_schur = f_p - C * Hinv_f_s;

    % Solve reduced system for primary variables
    x_p = S \ rhs_schur;

    % Reconstruct secondary variables
    x_s = Hinv_f_s - Hinv_B * x_p;

    % Assemble full solution
    sol = zeros(length(rhs), 1);
    sol(secondary_idx) = x_s;
    sol(primary_idx)   = x_p;
end

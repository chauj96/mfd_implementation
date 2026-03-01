function [m_proj, p_proj] = projectAnalyticalField3D(cell_struct, face_struct, phys, a, b, c, d)
% p(x,y,z) = a x + b y + c z + d
% gradp = [a; b; c]
%
% m_proj(f) = - A_f * (K * gradp) · n_f

    nCells = length(cell_struct);
    nFaces = length(face_struct);

    K_tensor = phys.K_tensor;     
    gradp = [a; b; c];             

    % Project pressure onto cell centers
    p_proj = zeros(nCells,1);
    for k = 1:nCells
        xc = cell_struct(k).center;        % [x y z]
        p_proj(k) = a*xc(1) + b*xc(2) + c*xc(3) + d;
    end

    % Project flux onto faces
    m_proj = zeros(nFaces,1);
    for f = 1:nFaces
        n_f = face_struct(f).normal(:);
        n_f = n_f / norm(n_f);
        A_f = face_struct(f).area;

        m_proj(f) = -A_f * dot(K_tensor * gradp, n_f);
    end
end

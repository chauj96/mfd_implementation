function [m_local_proj, p_local_proj, d_local_proj] = projectLocalAnalyticalField3D(cell_idx, cell_struct, face_struct, face_centers, a, b, c, d)
% p(x,y,z) = a x + b y + c z + d
% gradp = [a; b; c]
%
% m_proj(f) = - A_f * (K * gradp) · n_f

    % K_tensor = phys.K_tensor;     
    gradp = [a; b; c];             

    % Project pressure onto cell centers
    xc = cell_struct(cell_idx).center;
    K_tensor = cell_struct(cell_idx).K;
    p_local_proj = a*xc(1) + b*xc(2) + c*xc(3) + d;
    face_ids = cell_struct(cell_idx).faces;
    d_local_proj = a*face_centers(face_ids,1) + b*face_centers(face_ids,2) + c*face_centers(face_ids,3) + d;

    n_faces = length(face_ids);
    m_local_proj = zeros(n_faces,1);
    for idx = 1:n_faces
        f = face_ids(idx);
        n_f = face_struct(f).normal(:);
        n_f = n_f / norm(n_f);
        A_f = face_struct(f).area;
        m_local_proj(idx) = -A_f * dot(K_tensor * gradp, n_f);
    end

end

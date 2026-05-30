function [m_proj, p_proj] = projectAnalyticalField3D(cell_struct, face_struct, phys, a, b, c, d)
% p(x,y,z) = a x + b y + c z + d
% gradp = [a; b; c]
%
% m_proj(f) = - A_f * (K * gradp) · n_f

    nCells = length(cell_struct);
    nFaces = length(face_struct);

    % K_tensor = phys.K_tensor;     
    gradp = [a; b; c];             

    % Project pressure onto cell centers
    p_proj = zeros(nCells,1);
    for k = 1:nCells
        xc = cell_struct(k).center;     
        p_proj(k) = a*xc(1) + b*xc(2) + c*xc(3) + d;
    end

    % Project flux onto faces
    % K_tensor = phys.K_tensor;
    % m_proj = zeros(nFaces,1);
    % for f = 1:nFaces
    %     n_f = face_struct(f).normal(:);
    %     n_f = n_f / norm(n_f);
    %     A_f = face_struct(f).area;
    % 
    %     m_proj(f) = -A_f * dot(K_tensor * gradp, n_f);
    % end
   
    m_proj = zeros(nFaces,1);

    % for f = 1:nFaces
    % 
    %     n_f = face_struct(f).normal(:);
    %     n_f = n_f / norm(n_f);
    % 
    %     A_f = face_struct(f).area;
    % 
    %     neigh_cells = face_struct(f).cells;
    % 
    %     m_sum = 0;
    %     n_contrib = 0;
    % 
    %     for j = 1:length(neigh_cells)
    % 
    %         c_id = neigh_cells(j);
    % 
    %         if c_id <= 0
    %             continue;
    %         end
    % 
    %         Kc = cell_struct(c_id).K;
    % 
    %         mf_local = -A_f * dot(Kc * gradp, n_f);
    % 
    %         m_sum = m_sum + mf_local;
    %         n_contrib = n_contrib + 1;
    %     end
    % 
    %     m_proj(f) = m_sum / n_contrib;
    % end

    % Harmonic average
    for f = 1:nFaces

        n_f = face_struct(f).normal(:);
        n_f = n_f / norm(n_f);
    
        A_f = face_struct(f).area;
    
        neigh_cells = face_struct(f).cells;
        valid_cells = neigh_cells(neigh_cells > 0);
    
        % Boundary face
        if length(valid_cells) == 1
    
            Kf = cell_struct(valid_cells(1)).K;
    
        % Interior face
        else
    
            cL = valid_cells(1);
            cR = valid_cells(2);
    
            KL = cell_struct(cL).K;
            KR = cell_struct(cR).K;
    
            xL = cell_struct(cL).center(:);
            xR = cell_struct(cR).center(:);
            xf = face_struct(f).center(:);
    
            dL = norm(xf - xL);
            dR = norm(xf - xR);
    
            % Component-wise weighted harmonic average
            kL = diag(KL);
            kR = diag(KR);
    
            kf = (dL + dR) ./ (dL ./ kL + dR ./ kR);
    
            Kf = diag(kf);
    
        end
    
        m_proj(f) = -A_f * dot(Kf * gradp, n_f);
    
    end
end

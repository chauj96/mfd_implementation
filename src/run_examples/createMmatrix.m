function [cell_struct] = createMmatrix(cell_struct, face_struct, ip_type)
    
    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa','general_parametric','simple','quasi_tpfa','bdvlm'})}
    end
    
    n_cells = length(cell_struct);
    n_faces = length(face_struct);
    dim = length(face_struct(1).center);
    
    for c = 1:n_cells
        face_ids = cell_struct(c).faces;
        cell_nf = length(face_ids);
        Cc = cell_struct(c).center(:);
        K = cell_struct(c).K;
        v = cell_struct(c).volume;
        signs = double(cell_struct(c).faces_orientation(:));
    
        C = zeros(cell_nf, dim);
        N = zeros(cell_nf, dim);
        a = zeros(cell_nf, 1);
    
        for k = 1:cell_nf
            f = face_ids(k);
            Cf = face_struct(f).center(:);
            Nf = face_struct(f).normal(:);
            Af = face_struct(f).area;
    
            df = Cf - Cc;
            signf = sign(df' / norm(df) * Nf);
            assert(signf == signs(k), 'Orientation mismatch in cell %d', c);
    
            C(k,:) = df';
            N(k,:) = Af * signf * Nf';
            a(k) = Af;
        end
    
        if strcmp(ip_type,'tpfa')
    
            td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
            invT = diag(1 ./ abs(td));
    
        elseif strcmp(ip_type,'simple')
    
            t = 6 * sum(diag(K)) / dim;
            Q = orth(N ./ a);
            U = eye(cell_nf) - Q * Q';
            di = diag(1 ./ a);
            invT_reg = (v / t) * (di * U * di);
            invT = (C * (K \ C')) / v + invT_reg;
    
        elseif strcmp(ip_type,'quasi_tpfa')
    
            W = N * K * N';
            Qn = orth(N);
            P = eye(cell_nf) - Qn * Qn';
            diW = diag(1 ./ diag(W));
            invT_reg = (v / 2) * (P * diW * P);
            invT = (C * (K \ C')) / v + invT_reg;
    
        elseif strcmp(ip_type,'general_parametric')
    
            W = N * K * N';
            Qn = orth(N);
            P = eye(cell_nf) - Qn * Qn';
            diW = diag(1 ./ diag(W));
            invT_reg = (v / cell_nf) * (P * diW * P);
            invT = (C * (K \ C')) / v + invT_reg;
    
        elseif strcmp(ip_type,'bdvlm')
    
            R = diag(a) * C;
            Nbd = (N ./ a) * K;
            M0 = R * ((R' * Nbd) \ R');
            PN = eye(cell_nf) - Nbd * ((Nbd' * Nbd) \ Nbd');
            invT = diag(1 ./ a) * (M0 + (1/cell_nf) * PN) * diag(1 ./ a);
    
        end
    
        cell_struct(c).M = invT;
    
    end

end
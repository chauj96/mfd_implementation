function M = buildMmatrix(cell_struct, face_struct, ip_type)
    % Build the M matrix (inner product for velocity DOFs)
    %
    % ip_type: 'tpfa' or 'general_parametric'
    % t: number of faces per cell (required for 'general_parametric')

    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa', 'general_parametric'})}
    end

    rows = [];
    cols = [];
    vals = [];

    n_cells = length(cell_struct);
    n_faces = length(face_struct);

    dim = length(face_struct(1).center);

    switch ip_type
        case 'tpfa'
            for c = 1:n_cells
                face_ids = cell_struct(c).faces;
                for j = 1:length(face_ids)
                    f = face_ids(j);
                    xC = cell_struct(c).center(:);
                    xF = face_struct(f).center(:);
                    nf = face_struct(f).normal(:); 
                    Af = face_struct(f).area;  
                    K = nf' * cell_struct(c).K * nf; % extend to permeability tensor
                    d = abs(dot(xF - xC, nf));
                    invT = d / (K * Af);

                    rows(end+1) = f;
                    cols(end+1) = f;
                    vals(end+1) = invT;
                end
            end

        case 'general_parametric'
            for c = 1:n_cells
                
                face_ids = cell_struct(c).faces;
                cell_nf = length(face_ids);
                Cc = cell_struct(c).center(:);
                K = cell_struct(c).K;
                v = cell_struct(c).volume;

                % Local geometry arrays
                C = zeros(cell_nf, dim); % face direction vectors
                N = zeros(cell_nf, dim); % face normals
                a = zeros(cell_nf, 1);   % face areas
                for j = 1:cell_nf
                    f = face_ids(j);
                    Cf = face_struct(f).center(:);
                    Nf = face_struct(f).normal(:);
                    Af = face_struct(f).area;  

                    d = Cf - Cc; 
                    df = d;
                    signf = sign( df' * Nf );
                    C(j,:) = df';
                    N(j,:) = Af * signf * (Nf');
                    a(j)   = Af;
                end
                
                % general_parametric
                W = N * K * N';
                Q = orth(bsxfun(@rdivide, N, a));
                P = eye(cell_nf) - Q * Q';
                di = diag(1 ./ diag(W));
                invT = (C * (K \ C'))./v + (v / cell_nf) * (P * di * P);                

                % simple
                % t = 6 * sum(diag(K)) / size(K,2);
                % Q  = orth(bsxfun(@rdivide, N, a));
                % U  = eye(length(a)) - Q*Q';
                % di = diag(1 ./ a);
                %invT  = (C * (K \ C'))./v + (v / t)*(di * U * di);
                
                % tpfa
                td = sum(C .* (N*K), 2) ./ sum(C.*C, 2);
                %invT = diag(1 ./ abs(td));

                for i = 1:cell_nf
                    fi = face_ids(i);
                    for j = 1:cell_nf
                        fj = face_ids(j);
                        rows(end+1) = fi;
                        cols(end+1) = fj;
                        vals(end+1) = invT(i,j);
                    end
                end
            end
    end

    M = sparse(rows, cols, vals, n_faces, n_faces);
end

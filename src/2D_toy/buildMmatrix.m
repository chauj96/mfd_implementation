function M = buildMmatrix(cell_struct, face_struct, ip_type, t)
    % Build the M matrix (inner product for velocity DOFs)
    %
    % ip_type: 'tpfa' or 'general_parametric'
    % t: number of faces per cell (required for 'general_parametric')

    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa', 'general_parametric'})}
        t (1,1) double {mustBeNonnegative} = 6
    end

    rows = [];
    cols = [];
    vals = [];

    n_cells = length(cell_struct);
    n_faces = length(face_struct);

    % Automatically detect dimension (2D or 3D)
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
                    K = cell_struct(c).K;
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
                    

                    df = Cf - Cc;

                    C(j,:) = df';
                    N(j,:) = Af * sign(dot(df, Nf)) * Nf';
                    a(j)   = Af;
                end

                kappa = N * K * N';      % (nf x nf)
                Q = orth(N);             % (nf x rank(N))
                P = eye(cell_nf) - Q * Q';     % (nf x nf)
                di = diag(1 ./ diag(kappa));  % (nf x nf)

                invT = (C * inv(K) * C') ./ v + (v / t) .* (P * di * P);

                % TPFA version
                % td = sum(C .* (N*K), 2) ./ sum(C.*C, 2);
                % invT = diag(1 ./ abs(td));

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

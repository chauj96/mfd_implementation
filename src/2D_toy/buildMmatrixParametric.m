function M = buildMmatrixParametric(cell_struct, face_struct, ip_type)
    % Build the M matrix (inner product for velocity DOFs)
    %
    % ip_type: 'tpfa', 'general_parametric', or 'simple'

    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa', 'general_parametric', 'simple'})}
    end

    rows = [];
    cols = [];
    vals = [];

    n_cells = length(cell_struct);
    n_faces = length(face_struct);

    dim = length(face_struct(1).center);

    for c = 1:n_cells
        face_ids = cell_struct(c).faces;
        cell_nf = length(face_ids);
        Cc = cell_struct(c).center(:);
        K = cell_struct(c).K;
        v = cell_struct(c).volume;
        signs = cell_struct(c).faces_orientation;

        % Local geometry arrays
        C = zeros(cell_nf, dim); % face direction vectors
        N = zeros(cell_nf, dim); % face normals
        a = zeros(cell_nf, 1);   % face areas

        for idx = 1:cell_nf
            f = face_ids(idx);
            Cf = face_struct(f).center(:);
            Nf = face_struct(f).normal(:);
            Af = face_struct(f).area;

            d = Cf - Cc; 
            df = d;
            signf = sign(df' * Nf);
            C(idx,:) = df';
            N(idx,:) = Af * signf * (Nf');
            a(idx)   = Af;
        end

        % Compute inverse inner product matrix (invT) depending on ip_type
        switch ip_type
            case 'general_parametric'
                W = N * K * N';
                Q = orth(N);
                P = eye(cell_nf) - Q * Q';
                di = diag(1 ./ diag(W));
                invT_reg = (v / 4) * (P * di * P);
                invT = (C * (K \ C'))./v + invT_reg;

                % test regularization term
                Cr = norm(invT_reg * N);
                assert(Cr < 1.0e-12);

            case 'simple'
                t = 6 * sum(diag(K))/dim;
                Q = orth(bsxfun(@rdivide, N, a));
                U = eye(cell_nf) - Q * Q';
                di = diag(1 ./ a);
                invT_reg = (v / t)*(di * U * di);
                invT = (C * (K \ C'))./v + invT_reg;

                % test regularization term
                Cr = norm(invT_reg * N);
                assert(Cr < 1.0e-12);
                
                % Qs = orth(bsxfun(@times, C, a));
                % Us = eye(cell_nf) - Qs * Qs';
                % di = diag(a);
                % T = (N * K * N')./v + (t/v)*(di * Us * di);
                %invT = inv(T);
                %norm(inv(T) - invT)

            case 'tpfa'
                td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
                invT = diag(1 ./ abs(td));
        end

        c
        % test consistency conditions
        Cm = norm(invT * N * K - C);
        assert(Cm < 1.0e-12);

        % test geometrical property
        vol_res = norm(v*eye(dim) - C' * N);
        assert(vol_res < 1.0e-12);

        % invT
        % face_ids
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

    M = sparse(rows, cols, vals, n_faces, n_faces);

 % full(M)
% imagesc(M);
% colorbar;
% title('M Plot');
end

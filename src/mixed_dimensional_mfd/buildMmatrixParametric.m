function M = buildMmatrixParametric(cell_struct, face_struct, ip_type)

    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa', 'general_parametric', 'simple'})}
    end

    n_cells = length(cell_struct);
    n_faces = length(face_struct);
    dim = length(face_struct(1).center);

    total_nnz = sum(arrayfun(@(c) length(c.faces)^2, cell_struct));
    rows = zeros(total_nnz, 1);
    cols = zeros(total_nnz, 1);
    vals = zeros(total_nnz, 1);
    idx = 0;

    for c = 1:n_cells
        face_ids = cell_struct(c).faces;
        cell_nf = length(face_ids);
        Cc = cell_struct(c).center(:);
        K = cell_struct(c).K;
        v = cell_struct(c).volume;
        signs = cell_struct(c).faces_orientation;

        % Local C, N, a
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
            a(k)   = Af;
        end

        % ---------- Build invT based on initial ip_type ----------
        % Also compute R (using the initial choice)
        if strcmp(ip_type,'tpfa')
            % TPFA invT
            td   = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
            invT = diag(1 ./ abs(td));

        elseif strcmp(ip_type,'simple')
            % SIMPLE from start
            t  = 6 * sum(diag(K)) / dim;
            Q  = orth(N ./ a);
            U  = eye(cell_nf) - Q * Q';
            di = diag(1 ./ a);
            invT_reg = (v / t) * (di * U * di);
            invT = (C * (K \ C')) / v + invT_reg;

        else % 'general_parametric'
            W  = N * K * N';
            Qn = orth(N);
            P  = eye(cell_nf) - Qn * Qn';
            diW = diag(1 ./ diag(W));
            invT_reg = (v / cell_nf) * (P * diW * P);
            invT = (C * (K \ C')) / v + invT_reg;

        end

        % geometric consistency
        % vol_res = norm(v * eye(dim) - C' * N);
        % assert(vol_res < 1e-12);

        % assemble with orientation correction
        for i = 1:cell_nf
            fi = face_ids(i); si = signs(i);
            for j = 1:cell_nf
                fj = face_ids(j); sj = signs(j);
                idx = idx + 1;
                rows(idx) = fi;
                cols(idx) = fj;
                vals(idx) = si * sj * invT(i,j);
            end
        end
    end

    M = sparse(rows, cols, vals, n_faces, n_faces);
end



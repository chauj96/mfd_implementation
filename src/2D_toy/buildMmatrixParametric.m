function [M, Rvals, scheme_labels, cell_struct] = buildMmatrixParametric(cell_struct, face_struct, ip_type, tol)

    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa', 'general_parametric', 'simple'})}
        tol (1,1) double = 1
    end

    n_cells = length(cell_struct);
    n_faces = length(face_struct);
    dim = length(face_struct(1).center);

    Rvals = zeros(n_cells,1);
    scheme_labels = repmat({'M'}, n_cells, 1);

    total_nnz = sum(arrayfun(@(c) length(c.faces)^2, cell_struct));
    rows = zeros(total_nnz, 1);
    cols = zeros(total_nnz, 1);
    vals = zeros(total_nnz, 1);
    idx = 0;

    % ==== for error analysis!!! ======
    % TPFA_count = n_cells;
    % TPFA_threshold = 9888;
    % count = 0;
    % ==== for error analysis!!! ======

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

            % R from TPFA attempt
            s_mnk = svd(invT * N * K, 'econ');
            s_c   = svd(C, 'econ');
            k_mnk = max(s_mnk) / min(s_mnk);
            k_c   = max(s_c) / min(s_c);
            Rvals(c) = k_mnk / k_c;
            % scheme_labels{c} = 'T';

            % Decide: keep TPFA (T) or switch to SIMPLE (M)
            % if abs(Rvals(c) - 1) <= tol
            if abs(invT * N * K - C) < tol
                scheme_labels{c} = 'T';

            else
                % fprintf('TPFA is not consistent, R = %.6f\n', Rvals(c));
                % switch to SIMPLE MFD
                t  = 6 * sum(diag(K)) / dim;
                Q  = orth(N ./ a);
                U  = eye(cell_nf) - Q * Q';
                di = diag(1 ./ a);
                invT_reg = (v / t) * (di * U * di);
                invT = (C * (K \ C')) / v + invT_reg;
                scheme_labels{c} = 'M';
                % count = count + 1;
                % if count > 112
                %   td   = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
                %   invT = diag(1 ./ abs(td));
                %   scheme_labels{c} = 'T';
                % end
            end

        elseif strcmp(ip_type,'simple')
            % SIMPLE from start
            t  = 6 * sum(diag(K)) / dim;
            Q  = orth(N ./ a);
            U  = eye(cell_nf) - Q * Q';
            di = diag(1 ./ a);
            invT_reg = (v / t) * (di * U * di);
            invT = (C * (K \ C')) / v + invT_reg;

            % R using the current invT
            s_mnk = svd(invT * N * K, 'econ');
            s_c   = svd(C, 'econ');
            Rvals(c) = (max(s_mnk)/min(s_mnk)) / (max(s_c)/min(s_c));
            scheme_labels{c} = 'M';

        else % 'general_parametric'
            W  = N * K * N';
            Qn = orth(N);
            P  = eye(cell_nf) - Qn * Qn';
            diW = diag(1 ./ diag(W));
            invT_reg = (v / cell_nf) * (P * diW * P);
            invT = (C * (K \ C')) / v + invT_reg;

            s_mnk = svd(invT * N * K, 'econ');
            s_c   = svd(C, 'econ');
            Rvals(c) = (max(s_mnk)/min(s_mnk)) / (max(s_c)/min(s_c));
            scheme_labels{c} = 'M';
        end

        % geometric consistency
        vol_res = norm(v * eye(dim) - C' * N);
        assert(vol_res < 1e-12);

        % store R value directly in cell struct
        cell_struct(c).R = Rvals(c);

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


% function M = buildMmatrixParametric(cell_struct, face_struct, ip_type)
%     % Build the M matrix (inner product for velocity DOFs)
%     % ip_type: 'tpfa', 'general_parametric', or 'simple'
% 
%     arguments
%         cell_struct
%         face_struct
%         ip_type (1,:) char {mustBeMember(ip_type, {'tpfa', 'general_parametric', 'simple'})}
%     end
% 
%     n_cells = length(cell_struct);
%     n_faces = length(face_struct);
%     dim = length(face_struct(1).center);
% 
%     % === Precompute exact number of non-zeros ===
%     total_nnz = sum(arrayfun(@(c) length(c.faces)^2, cell_struct));
%     rows = zeros(total_nnz, 1);
%     cols = zeros(total_nnz, 1);
%     vals = zeros(total_nnz, 1);
%     idx = 0;
% 
%     for c = 1:n_cells
%         face_ids = cell_struct(c).faces;
%         cell_nf = length(face_ids);
%         Cc = cell_struct(c).center(:);
%         K = cell_struct(c).K;
%         v = cell_struct(c).volume;
%         signs = cell_struct(c).faces_orientation;
% 
%         % Local arrays
%         C = zeros(cell_nf, dim); 
%         N = zeros(cell_nf, dim); 
%         a = zeros(cell_nf, 1);   
% 
%         for k = 1:cell_nf
%             f = face_ids(k);
%             Cf = face_struct(f).center(:);
%             Nf = face_struct(f).normal(:);
%             Af = face_struct(f).area;
% 
%             df = Cf - Cc;
%             signf = sign(df' / norm(df) * Nf);
%             assert(signf == signs(k), 'Orientation mismatch in cell %d', c);
% 
%             C(k,:) = df';
%             N(k,:) = Af * signf * Nf';
%             a(k)   = Af;
%         end
% 
%         % Compute inner product
%         switch ip_type
%             case 'general_parametric'
%                 W = N * K * N';
%                 Q = orth(N);
%                 P = eye(cell_nf) - Q * Q';
%                 di = diag(1 ./ diag(W));
%                 invT_reg = (v / cell_nf) * (P * di * P);
%                 invT = (C * (K \ C')) / v + invT_reg;
% 
%                 s_mnk = svd(invT * N * K, 'econ');
%                 s_c = svd(C, 'econ');
%                 k_mnk = max(s_mnk) / min(s_mnk);
%                 k_c = max(s_c) / min(s_c);
% 
%                 R = k_mnk / k_c
% 
%                 assert(norm(invT_reg * N) < 1e-12);
%                 assert(norm(invT * N * K - C) < 1e-12);
% 
%             case 'simple'
%                 t = 6 * sum(diag(K)) / dim;
%                 Q = orth(N ./ a);
%                 U = eye(cell_nf) - Q * Q';
%                 di = diag(1 ./ a);
%                 invT_reg = (v / t) * (di * U * di);
%                 invT = (C * (K \ C')) / v + invT_reg;
% 
%                 s_mnk = svd(invT * N * K, 'econ');
%                 s_c = svd(C, 'econ');
%                 k_mnk = max(s_mnk) / min(s_mnk);
%                 k_c = max(s_c) / min(s_c);
% 
%                 R = k_mnk / k_c;
% 
%                 assert(norm(invT_reg * N) < 1e-12);
%                 assert(norm(invT * N * K - C) < 1e-12);
% 
%             case 'tpfa'
%                 td = sum(C .* (N * K), 2) ./ sum(C .* C, 2);
%                 invT = diag(1 ./ abs(td));
%                 s_mnk = svd(invT * N * K, 'econ');
%                 s_c = svd(C, 'econ');
% 
%                 k_mnk = max(s_mnk) / min(s_mnk);
%                 k_c = max(s_c) / min(s_c);
% 
%                 R = k_mnk / k_c;
% 
%                 % if norm(invT * N * K - C) > 1e-12
%                 if norm(R - 1) > 1e-12
%                     warning('TPFA not consistent at cell %d. Switch to Simple', c);
%                     t = 6 * sum(diag(K)) / dim;
%                     Q = orth(N ./ a);
%                     U = eye(cell_nf) - Q * Q';
%                     di = diag(1 ./ a);
%                     invT_reg = (v / t) * (di * U * di);
%                     invT = (C * (K \ C')) / v + invT_reg;
%                 end
%         end
% 
%         % Geometric consistency check
%         vol_res = norm(v * eye(dim) - C' * N);
%         assert(vol_res < 1e-12);
% 
%         % Assemble with orientation correction
%         for i = 1:cell_nf
%             fi = face_ids(i);
%             si = signs(i);
%             for j = 1:cell_nf
%                 fj = face_ids(j);
%                 sj = signs(j);
% 
%                 idx = idx + 1;
%                 rows(idx) = fi;
%                 cols(idx) = fj;
%                 vals(idx) = si * sj * invT(i,j);
%             end
%         end
%     end
% 
%     M = sparse(rows, cols, vals, n_faces, n_faces);
% end

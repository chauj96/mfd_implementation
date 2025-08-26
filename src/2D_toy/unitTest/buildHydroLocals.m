function locals = buildHydroLocals(cell_struct, face_struct, ip_type)
    arguments
        cell_struct
        face_struct
        ip_type (1,:) char {mustBeMember(ip_type, {'tpfa','general_parametric','simple'})}
    end

    n_cells = length(cell_struct);
    dim = length(face_struct(1).center);
    locals(n_cells) = struct('faces',[],'Cloc',[],'Nloc',[],'Mloc',[],'center',[]);

    for c = 1:n_cells
        face_ids = cell_struct(c).faces;
        nf = length(face_ids);
        Cc = cell_struct(c).center(:);
        K = cell_struct(c).K;
        V = cell_struct(c).volume;
        signs = cell_struct(c).faces_orientation;

        C = zeros(nf, dim); 
        N = zeros(nf, dim); 
        a = zeros(nf, 1);

        for k = 1:nf
            f  = face_ids(k);
            Cf = face_struct(f).center(:);
            Nf = face_struct(f).normal(:);
            Af = face_struct(f).area;

            df = Cf - Cc;
            signf = sign((df'/norm(df)) * Nf);
            assert(signf == signs(k), 'Orientation mismatch in cell %d', c);

            C(k,:) = df';
            N(k,:) = (Af * signf) * Nf';
            a(k) = Af;
        end

        % === build transmissibility matrix T (not M for this time) ===
        switch ip_type
            case 'general_parametric'
                t = 2;    % QTPFA when t = 2
                W = N * K * N.';               
                QC = orth(C);                     
                PC = eye(size(C, 1)) - QC * QC.';   
                Tloc = (W + t * PC * diag(diag(W)) * PC) / V;

            case 'simple'
                A  = diag(a);
                W = N * K * N';     
                AC = A * C;           
                S  = C.' * (A * A) * C;  
                U  = eye(nf) - AC * (S \ AC');   
                Treg = (6 / dim) * trace(K) * A * U * A;
                Tloc = (W + Treg) / V;

            case 'tpfa'
                td   = sum(C .* (N*K), 2) ./ sum(C .* C, 2);
                Tloc = diag(abs(td));

        end

        % geometric identity check
        assert(norm(V*eye(dim) - C'*N, 'fro') < 1e-12);

        locals(c).faces  = face_ids;
        locals(c).Cloc   = C;
        locals(c).Nloc   = N;   
        locals(c).Tloc   = Tloc;
        locals(c).center = Cc;
    end
end

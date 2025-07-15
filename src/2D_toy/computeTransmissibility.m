function face_struct = computeTransmissibility(cell_struct, face_struct)
% Compute transmissibility for each face based on harmonic averaging

    for f = 1:length(face_struct)
        cids = face_struct(f).cells;       % cell indices (could be 1 or 2)
        nf = face_struct(f).normal(:);     % face normal (column vector)
        Af = face_struct(f).area;          % face area
    
        if length(cids) == 2
            % Interior face
            cL = cids(1);
            cR = cids(2);
    
            xL = cell_struct(cL).center(:);
            xR = cell_struct(cR).center(:);
    
            % permeability
            K1 = cell_struct(cL).K;
            K2 = cell_struct(cR).K;
    
            % normal projection
            d = abs(dot(xR - xL, nf));
    
            % harmonic average + distance corrector
            K_harm = 2 * K1 * K2 / (K1 + K2);
            T = K_harm * Af;
    
        else
            % Boundary face
            c = cids(1);
            xC = cell_struct(c).center(:);
            xF = face_struct(f).center(:);
    
            K = cell_struct(c).K;
    
            d = abs(dot(xF - xC, nf));
            T = K * Af;
        end
    
        face_struct(f).T = T;
    end


end
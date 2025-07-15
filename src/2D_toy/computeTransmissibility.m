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
            xF = face_struct(f).center(:);

            % permeability
            KL = cell_struct(cL).K;
            KR = cell_struct(cR).K;
            
    
            % normal projection
            dL = abs(dot(xF - xL, nf));
            dR = abs(dot(xF - xR, nf));
            invTL = dL / KL;
            invTR = dR / KR;
            invT = Af * (invTL+invTR);
    
        else
            % Boundary face
            c = cids(1);
            xC = cell_struct(c).center(:);
            xF = face_struct(f).center(:);
    
            K = cell_struct(c).K;
    
            d = abs(dot(xF - xC, nf));
            invT = d / (K * Af);
        end
    
        face_struct(f).invT = invT;
    end


end
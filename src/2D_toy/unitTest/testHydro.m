function testHydro(cell_struct, face_struct, locals, rho, gvec, pD, tol)
    if nargin < 7
        tol = 1e-10; 
    end
    
    isB = arrayfun(@(f) length(f.cells) == 1, face_struct);
    [~, I] = max(arrayfun(@(f) f.center(2), face_struct));
    xTop = face_struct(I).center(:);
    p0 = pD + rho * (gvec.' * xTop);
    
    maxErr = 0; 
    relErr = 0;
    for c = 1:numel(cell_struct)
        xE = locals(c).center(:);
        fids = locals(c).faces; 
        nf = numel(fids);
    
        pE = p0 - rho * (gvec.' * xE);
        piE = zeros(nf, 1);
        for k = 1:nf
            xf = face_struct(fids(k)).center(:);
            piE(k) = p0 - rho * (gvec.' * xf);
        end
    
        CE = locals(c).Cloc;
        bE = CE * (rho * gvec);    
    
        % u = T(pi - e p) + T C(rho g)
        Tloc = locals(c).Tloc;
        uE = Tloc * (piE - pE * ones(nf, 1) + bE);
    
        maxErr = max(maxErr, max(abs(uE)));
        relErr  = max(relErr, norm(uE) / max(1,norm(bE)));
    end
    
fprintf('maxErr = %.3e, relErr = %.3e\n', maxErr, relErr);
assert(maxErr <= tol, 'Hydro well-balanced (u=0) failed');
end





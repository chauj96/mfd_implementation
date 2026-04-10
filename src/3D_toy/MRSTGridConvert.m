function [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G)
% Convert MRST 3D grid into custom format compatible with existing MFD code.
%
% NOTE: G must already have geometry computed (via computeGeometry) and
% G.cells.centroids must contain the VORONOI GENERATORS (sites), not the
% geometric centroids.  Do NOT call computeGeometry inside this function
% as it would overwrite the generators with geometric centroids, breaking
% the MFD consistency condition C'*N = v*K.
%
% Both scratch() and scratch_conformal() guarantee this on exit.

    V3 = G.nodes.coords;
    nCells = G.cells.num;
    nFaces = G.faces.num;

    cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'face_normals', {}, 'volume', {});
    face_struct = struct('cells', {}, 'verts', {}, 'center', {}, 'normal', {}, 'area', {});
    cells3D = cell(nCells,1);

    % Build face_struct
    faceNodePos = G.faces.nodePos;
    faceNodes   = G.faces.nodes;

    for f = 1:nFaces
        ids   = faceNodePos(f):faceNodePos(f+1)-1;
        verts = faceNodes(ids);

        neigh = G.faces.neighbors(f,:);
        neigh = neigh(neigh > 0);

        n_unit = G.faces.normals(f,:)'/G.faces.areas(f);

        face_struct(f).cells  = neigh(:)';
        face_struct(f).verts  = verts(:);
        face_struct(f).center = G.faces.centroids(f,:)';
        face_struct(f).normal = n_unit(:);
        face_struct(f).area   = G.faces.areas(f);
    end

    % Build cell_struct
    cf   = G.cells.faces(:,1);
    cpos = G.cells.facePos;

    for c = 1:nCells
        ids        = cpos(c):cpos(c+1)-1;
        cell_faces = cf(ids);

        % Use Voronoi generator as stored center (needed for MFD consistency)
        cell_struct(c).center = G.cells.centroids(c,:)';
        cell_struct(c).volume = G.cells.volumes(c);
        cell_struct(c).faces  = cell_faces(:)';

        % Compute a GEOMETRIC reference point (mean of face centroids) for
        % orientation sign detection.  The Voronoi generator (cell center)
        % may lie outside the cell for boundary cells, making dot-product
        % sign tests unreliable.
        face_centers = zeros(numel(cell_faces), 3);
        for k = 1:numel(cell_faces)
            face_centers(k,:) = face_struct(cell_faces(k)).center';
        end
        xc_geom = cell_struct(c).center;   % geometric interior reference

        face_signs   = zeros(numel(cell_faces), 1);
        face_normals = zeros(numel(cell_faces), 3);
        verts_local  = [];

        for k = 1:numel(cell_faces)
            f  = cell_faces(k);

            nf = face_struct(f).normal(:);
            xf = face_struct(f).center(:);

            d = dot(nf, xf - xc_geom);   % use geometric center, not Voronoi site
            if abs(d) < 1e-14 * norm(xf - xc_geom) * norm(nf)
                warning('Degenerate face orientation for cell %d, face %d', c, f);
                s        = 1;
                nf_local = nf;
            elseif d < 0
                s        = -1;
                nf_local = -nf;
            else
                s        =  1;
                nf_local =  nf;
            end

            face_signs(k)      = s;
            face_normals(k,:)  = nf_local';

            verts_local = [verts_local; face_struct(f).verts(:)];
        end

        cell_struct(c).faces_orientation = face_signs(:)';
        cell_struct(c).face_normals      = face_normals;
        cells3D{c} = unique(verts_local(:), 'stable')';
    end

    fprintf('Converted MRST 3D grid:\n');
    fprintf('  %d vertices\n', size(V3,1));
    fprintf('  %d cells\n',    nCells);
    fprintf('  %d faces\n',    nFaces);

    % Run volume consistency check
    checkVolumeConsistency(cell_struct, face_struct, nCells);
end


function checkVolumeConsistency(cell_struct, face_struct, nCells)
% checkVolumeConsistency  Verify cell geometry via the MFD consistency condition.
%
%   Computes the tensor T = sum_f [ (s * n_f * A_f) * (x_f - x_c)^T ].
%   For a valid polyhedron, T must equal V_c * I (where I is 3x3 identity).

    inexact     = 0;
    max_rel_err = 0;
    fail_list   = [];

    fprintf('\n=== MRSTGridConvert:: MFD Volume & Tensor Consistency Check ===\n');

    for c = 1:nCells
        cell_faces = cell_struct(c).faces;
        face_signs = cell_struct(c).faces_orientation;
        Vc         = cell_struct(c).volume;
        xc         = cell_struct(c).center(:);

        T = zeros(3,3);

        for k = 1:numel(cell_faces)
            f  = cell_faces(k);
            s  = face_signs(k); % +1 or -1 outward sign
            xf = face_struct(f).center(:);
            nf = face_struct(f).normal(:); % unit normal
            af = face_struct(f).area;

            n_out = s * nf; % Outward area-weighted normal (3x1)
            T     = T + af * n_out * (xf - xc)';
        end

        % Trace/div-theorem volume
        V_div = trace(T) / 3;

        % Check MFD tensor equation
        rel_err     = norm(T - Vc * eye(3), 'fro') / max(abs(Vc), eps);
        rel_vol_err = abs(V_div - Vc) / max(abs(Vc), eps);
        max_rel_err = max(max_rel_err, rel_err);

        if rel_err > 1e-10
            inexact = inexact + 1;
            fail_list(end+1,:) = [c, rel_err, Vc, V_div, rel_vol_err]; 
        end
    end

    fprintf('  Cells checked           : %d\n', nCells);
    fprintf('  Max relative tensor err : %.3e\n', max_rel_err);

    if inexact == 0
        fprintf('  All cells PASS (tol 1e-10).\n');
    else
        fprintf('  FAILED cells (rel_err > 1e-10): %d / %d\n', inexact, nCells);

        vol_err = fail_list(:,5);
        fprintf('\n  Volume discrepancy (|V_div - V_stored|/V_stored) in failed cells:\n');
        fprintf('    min/mean/max: %.3e / %.3e / %.3e\n', ...
                min(vol_err), mean(vol_err), max(vol_err));

        [~, ord] = sort(fail_list(:,2), 'descend');
        fprintf('\n  Worst 5 cells (ci, rel_T_err, V_stored, V_div):\n');
        for k = 1:min(5, size(fail_list,1))
            r = fail_list(ord(k),:);
            fprintf('    cell %4d  rel_err=%.3e  V_stored=%.4e  V_div=%.4e\n', ...
                    r(1), r(2), r(3), r(4));
        end
    end
    fprintf('-----------------------------------------------------\n\n');
end

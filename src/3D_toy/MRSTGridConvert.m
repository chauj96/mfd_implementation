function [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G)
% We convert MRST 3D grid into custom format compatible with existing code

    G = computeGeometry(G);

    V3 = G.nodes.coords; 
    nCells = G.cells.num;
    nFaces = G.faces.num;

    cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'face_normals', {}, 'volume', {});
    face_struct = struct('cells', {}, 'verts', {}, 'center', {}, 'normal', {}, 'area', {});
    cells3D = cell(nCells,1);

    % Build face_struct
    faceNodePos = G.faces.nodePos;
    faceNodes = G.faces.nodes;

    for f = 1:nFaces
        ids = faceNodePos(f):faceNodePos(f+1)-1;
        verts = faceNodes(ids);

        neigh = G.faces.neighbors(f,:);
        neigh = neigh(neigh > 0);

        nvec = G.faces.normals(f,:)';
        nrm = norm(nvec);
        if nrm > 0
            n_unit = nvec / nrm;
        else
            n_unit = nvec;
        end

        face_struct(f).cells = neigh(:)';
        face_struct(f).verts = verts(:);
        face_struct(f).center = G.faces.centroids(f,:)';
        face_struct(f).normal = n_unit(:);
        face_struct(f).area = G.faces.areas(f);
    end

    % Build cell_struct
    cf = G.cells.faces(:,1);
    cpos = G.cells.facePos;

    for c = 1:nCells
        ids = cpos(c):cpos(c+1)-1;
        cell_faces = cf(ids);

        cell_struct(c).center = G.cells.centroids(c,:)';
        cell_struct(c).volume = G.cells.volumes(c);
        cell_struct(c).faces = cell_faces(:)';

        face_signs = zeros(numel(cell_faces),1);
        face_normals = zeros(numel(cell_faces),3); 
        verts_local = [];

        for k = 1:numel(cell_faces)
            f = cell_faces(k);

            nf = face_struct(f).normal(:);
            xf = face_struct(f).center(:);
            xc = cell_struct(c).center(:);

            if dot(nf, xf - xc) < 0
                nf_local = -nf;
                s = -1;
            else
                nf_local = nf;
                s = 1;
            end

            face_signs(k) = s;
            face_normals(k,:) = nf_local'; 

            verts_local = [verts_local; face_struct(f).verts(:)];
        end

        cell_struct(c).faces_orientation = face_signs(:)';
        cell_struct(c).face_normals = face_normals; 
        cells3D{c} = unique(verts_local, 'stable')';
    end

    fprintf('Converted MRST 3D grid:\n');
    fprintf('  %d vertices\n', size(V3,1));
    fprintf('  %d cells\n', nCells);
    fprintf('  %d faces\n', nFaces);
end

% function [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G)
% % We convert MRST 3D grid into custom format compatible with existing code
% 
%     G = computeGeometry(G);
% 
%     V3 = G.nodes.coords; 
%     nCells = G.cells.num;
%     nFaces = G.faces.num;
% 
%     cell_struct = struct('center', {}, 'faces', {}, 'faces_orientation', {}, 'volume', {});
%     face_struct = struct('cells', {}, 'verts', {}, 'center', {}, 'normal', {}, 'area', {});
%     cells3D = cell(nCells,1);
% 
%     % Build face_struct
%     faceNodePos = G.faces.nodePos;
%     faceNodes = G.faces.nodes;
% 
%     for f = 1:nFaces
%         ids = faceNodePos(f):faceNodePos(f+1)-1;
%         verts = faceNodes(ids);
% 
%         neigh = G.faces.neighbors(f,:);
%         neigh = neigh(neigh > 0);
% 
%         nvec = G.faces.normals(f,:)';
%         nrm = norm(nvec);
%         if nrm > 0
%             n_unit = nvec / nrm;
%         else
%             n_unit = nvec;
%         end
% 
%         face_struct(f).cells = neigh(:)';
%         face_struct(f).verts = verts(:);
%         face_struct(f).center = G.faces.centroids(f,:)';
%         face_struct(f).normal = n_unit(:);
%         face_struct(f).area = G.faces.areas(f);
%     end
% 
%     % Build cell_struct
%     cf = G.cells.faces(:,1);
%     cpos = G.cells.facePos;
% 
%     for c = 1:nCells
%         ids = cpos(c):cpos(c+1)-1;
%         cell_faces = cf(ids);
% 
%         cell_struct(c).center = G.cells.centroids(c,:)';
%         cell_struct(c).volume = G.cells.volumes(c);
%         cell_struct(c).faces = cell_faces(:)';
% 
%         face_signs = zeros(numel(cell_faces),1);
%         verts_local = [];
% 
%         for k = 1:numel(cell_faces)
%             f = cell_faces(k);
% 
%             nf = face_struct(f).normal(:);
%             xf = face_struct(f).center(:);
%             xc = cell_struct(c).center(:);
% 
%             s = sign(dot(nf, xf - xc));
%             if s == 0
%                 s = 1;
%             end
%             face_signs(k) = s;
% 
%             verts_local = [verts_local; face_struct(f).verts(:)];
%         end
% 
%         cell_struct(c).faces_orientation = face_signs(:)';
%         cells3D{c} = unique(verts_local, 'stable')';
%     end
% 
%     fprintf('Converted MRST 3D grid:\n');
%     fprintf('  %d vertices\n', size(V3,1));
%     fprintf('  %d cells\n', nCells);
%     fprintf('  %d faces\n', nFaces);
% end
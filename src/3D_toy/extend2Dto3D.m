function [cell3D, face3D, V3, cells3D] = extend2Dto3D(V2, cells2D, H)

    Nv = size(V2,1);
    
    % 3D vertices
    V0 = [V2, zeros(Nv,1)];
    V1 = [V2, H*ones(Nv,1)];
    V3 = [V0; V1];
    
    nCells = numel(cells2D);
    
    % build prism cell connectivity (bot/top)
    cells3D = cell(nCells,1);
    for c = 1:nCells
        bot = cells2D{c}(:)';
        top = bot + Nv;
        cells3D{c} = [bot, top];
    end
    
    % Global face map: key
    face_map = containers.Map('KeyType','char','ValueType','int32');
    face_verts = {};           
    face_cells = {};              
    
    cell_faces = cell(nCells,1);  % store face ids per cell
    
    % helper to register face (order-insensitive unique)
        function fid = registerFace(vlist, cell_id)
            v_sorted = sort(vlist);
            key = sprintf('%d_', v_sorted);
            if isKey(face_map, key)
                fid = face_map(key);
                face_cells{fid}(end+1) = cell_id;
            else
                fid = numel(face_verts) + 1;
                face_map(key) = fid;
                face_verts{fid} = vlist(:)';
                face_cells{fid} = cell_id;
            end
        end
    
    % Build faces per cell (bottom, top, sides)
    for c = 1:nCells
        bot = cells2D{c}(:)';
        top = bot + Nv;
        n = numel(bot);
    
        fids = zeros(n+2,1);
    
        % bottom face (polygon)
        fids(1) = registerFace(bot, c);
    
        % top face (polygon)
        fids(2) = registerFace(top, c);
    
        % side faces: each edge -> quad [vi, vj, vj+Nv, vi+Nv]
        for k = 1:n
            i = bot(k);
            j = bot(mod(k,n)+1);
            quad = [i, j, j+Nv, i+Nv];
            fids(2+k) = registerFace(quad, c);
        end
    
        cell_faces{c} = fids;
    end
    
    % Build face_struct geometry (center, normal, area)
    nFaces = numel(face_verts);
    face3D = struct('verts',{},'center',{},'normal',{},'area',{},'cells',{});
    
    for f = 1:nFaces
        vlist = face_verts{f};
        pts = V3(vlist,:);       
        fc = mean(pts,1);
    
        tol_axis = 1e-12;  
    
        xspan = max(pts(:,1)) - min(pts(:,1));
        yspan = max(pts(:,2)) - min(pts(:,2));
        zspan = max(pts(:,3)) - min(pts(:,3));
        
        if zspan < tol_axis
            n_unit = [0 0 1];     % bottom/top -> +ez
        elseif xspan < tol_axis
            n_unit = [1 0 0];     % side faces normal to x -> +ex
        elseif yspan < tol_axis
            n_unit = [0 1 0];     % side faces normal to y -> +ey
        else
            % fallback
            nsum = [0 0 0];
            for k = 2:(size(pts,1)-1)
                nsum = nsum + cross(pts(k,:) - pts(1,:), pts(k+1,:) - pts(1,:));
            end
            n_unit = nsum / norm(nsum);
        end
    
        % area = sum triangle areas
        area = 0;
        for k = 2:(size(pts,1)-1)
            area = area + 0.5 * norm(cross(pts(k,:) - pts(1,:), pts(k+1,:) - pts(1,:)));
        end
    
        face3D(f).verts = vlist(:);
        face3D(f).center = fc(:);
        face3D(f).normal = n_unit(:); 
        face3D(f).area = area;
        face3D(f).cells = face_cells{f}(:);
    end
    
    % Build cell_struct (center, volume, orientation signs)
    cell3D = struct('center',{},'faces',{},'faces_orientation',{},'volume',{});
    
    for c = 1:nCells
        bot = cells2D{c}(:);
        pts2 = V2(bot,:);
    
        % 2D polygon area (shoelace) -> volume = area*H
        x = pts2(:,1); y = pts2(:,2);
        A = 0.5*abs(sum(x.*circshift(y,-1)) - sum(y.*circshift(x,-1)));
        vol = A * H;
    
        % 3D cell center = average of its 3D vertices
        vids3 = cells3D{c};
        cc = mean(V3(vids3,:), 1);
    
        fids = cell_faces{c};
        sgns = zeros(numel(fids),1);
    
        for k = 1:numel(fids)
            fid = fids(k);
            nf = face3D(fid).normal(:);
            fc = face3D(fid).center(:);
            s = sign(dot(nf, fc - cc(:)));
            if s == 0, s = 1; end
            sgns(k) = s;
        end
    
        cell3D(c).center = cc(:);
        cell3D(c).faces = fids(:);
        cell3D(c).faces_orientation = sgns(:);
        cell3D(c).volume = vol;
    end

end

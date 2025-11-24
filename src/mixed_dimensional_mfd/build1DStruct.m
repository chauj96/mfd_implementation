function cell_struct_1d = build1DStruct(vertices, cell_struct_2d, nz, Lx, Lz)
% Return fracture information: center, length, faces, faces orientation,
% adjacent cells (parents cells) and permeability tensor (ortho,
% tangential)

    dz = Lz / nz;

    % Fracture is in the middle of the domain 
    fracture_x = Lx/2;

    idx = 0;

    for k = 1:nz

        idx = idx + 1;

        % Segment center
        yc = (k - 0.5) * dz;
        center = [fracture_x; yc];
        
        % Segment endpoints
        y1 = (k-1)*dz;
        y2 = k*dz;

        % Find vertex indices
        v_lower = find(abs(vertices(:,1)-fracture_x)<1e-12 & abs(vertices(:,2)-y1)<1e-12);
        v_upper = find(abs(vertices(:,1)-fracture_x)<1e-12 & abs(vertices(:,2)-y2)<1e-12);

        verts = [v_lower, v_upper];

        % Compute length
        len = dz;

        % -------- parent cells --------
        xs = cellfun(@(c)c(1), {cell_struct_2d.center});
        zs = cellfun(@(c)c(2), {cell_struct_2d.center});

        % left
        left_candidates = find(xs < fracture_x);
        % [~,iL] = min(abs(zs(left_candidates) - fracture_x));
        diffz = abs(zs(left_candidates) - yc);
        iL = find(diffz == min(diffz), 1, 'last');
        left_cell = left_candidates(iL);


        % right
        right_candidates = find(xs > fracture_x);
        [~,iR] = min(abs(zs(right_candidates) - yc));
        right_cell = right_candidates(iR);

        % Populate
        cell_struct_1d(idx).center = center;
        cell_struct_1d(idx).volume = len;
        cell_struct_1d(idx).faces = verts; % [vertex_bottom, vertex_top]
        cell_struct_1d(idx).faces_orientation = [-1, 1]; % [bottom, top]
        cell_struct_1d(idx).parent_2d = [left_cell, right_cell]; 

        % tangential / normal permeability
        K_ortho = 1.0; 
        K_tangent = 1.0; 
        cell_struct_1d(idx).K = [K_ortho, 0; 0, K_tangent];

    end
end


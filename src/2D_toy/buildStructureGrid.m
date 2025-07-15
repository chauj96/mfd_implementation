function [cell_struct, face_struct] = buildStructureGrid(nx, nz, Lx, Lz)
    
    dx = Lx / nx;
    dz = Lz / nz;
    
    n_cells = nx * nz;
    n_faces_x = (nx + 1) * nz;
    n_faces_z = (nz + 1) * nx;
    n_faces = n_faces_x + n_faces_z;
    
    % Face storage
    face_struct = struct('cells', {}, 'normal', {}, 'center', {}, 'area', {});
    face_counter = 0;
    face_map = containers.Map;
    
    % Vertical faces (left-right)
    for j = 1:nz
        for i = 1:(nx+1)
            face_counter = face_counter + 1;
            x = (i - 1) * dx;
            z = (j - 1) * dz + dz / 2;
    
            face_struct(face_counter).center = [x, z];
            face_struct(face_counter).normal = [1; 0];
            face_struct(face_counter).cells = [];
            face_struct(face_counter).area = dz;
    
            key = sprintf('v_%d_%d', i, j);
            face_map(key) = face_counter;
        end
    end
    
    % Horizontal faces (bottom-top)
    for j = 1:(nz+1)
        for i = 1:nx
            face_counter = face_counter + 1;
            x = (i - 1) * dx + dx / 2;
            z = (j - 1) * dz;
    
            face_struct(face_counter).center = [x, z];
            face_struct(face_counter).normal = [0; 1];
            face_struct(face_counter).cells = [];
            face_struct(face_counter).area = dx;
    
            key = sprintf('h_%d_%d', i, j);
            face_map(key) = face_counter;
        end
    end
    
    % Cells and their connection to faces
    cell_struct = struct('center', {}, 'faces', {}, 'face_dirs', {}, 'volume', {});
    cell_id = 0;
    
    for j = 1:nz
        for i = 1:nx
            cell_id = cell_id + 1;
    
            xc = (i - 1) * dx + dx / 2;
            zc = (j - 1) * dz + dz / 2;
            cell_struct(cell_id).center = [xc, zc];
            cell_struct(cell_id).volume = dx * dz;
    
            % Get face indices
            fL = face_map(sprintf('v_%d_%d', i, j));
            fR = face_map(sprintf('v_%d_%d', i+1, j));
            fB = face_map(sprintf('h_%d_%d', i, j));
            fT = face_map(sprintf('h_%d_%d', i, j+1));
    
            cell_struct(cell_id).faces = [fL, fR, fB, fT];
            cell_struct(cell_id).face_dirs = [-1, 1, -1, 1];
    
            % Update each face's connected cells
            face_struct(fL).cells = [face_struct(fL).cells, cell_id];
            face_struct(fR).cells = [face_struct(fR).cells, cell_id];
            face_struct(fB).cells = [face_struct(fB).cells, cell_id];
            face_struct(fT).cells = [face_struct(fT).cells, cell_id];
        end
    end

end

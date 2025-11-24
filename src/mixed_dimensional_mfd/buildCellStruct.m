function [cell_struct_2d, face_struct_2d] = buildCellStruct(cell_struct_2d, face_struct_2d, Lx, nx)

    fracture = Lx/2;
    dx = Lx/nx;
    n_cells = length(cell_struct_2d);

    all_faces = [cell_struct_2d.faces]; 
    max_face_id = max(all_faces);

    for c = 1:n_cells

        xc = cell_struct_2d(c).center(1);

        if abs(xc - fracture) > dx
            continue;
        end

        % Determine if the face is placed left or right of fracture
        is_right = xc > fracture;

        % face ordering: bottom(1), right(2), top(3), left(4)
        if ~is_right
            face_idx = 2;    
        else
            face_idx = 4;   
        end

        old_face = cell_struct_2d(c).faces(face_idx);

        if ~is_right
            % =======================
            % CASE 1: LEFT SIDE CELL
            % =======================
            face_struct_2d(old_face).cells = c;
        
        else
            % =======================
            % CASE 2: RIGHT SIDE CELL
            % =======================

            new_face = max_face_id + 1;
            max_face_id = new_face;

            cell_struct_2d(c).faces(face_idx) = new_face;

            face_struct_2d(new_face).center = face_struct_2d(old_face).center;
            face_struct_2d(new_face).normal = face_struct_2d(old_face).normal;
            face_struct_2d(new_face).area = face_struct_2d(old_face).area;
            face_struct_2d(new_face).cells = c;

        end
    end
end


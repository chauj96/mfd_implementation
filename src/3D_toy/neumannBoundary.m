function [A, rhs_BC, f_ids] = neumannBoundary(A, ~, face_struct)

    n_faces = length(face_struct);
    rhs_BC = zeros(n_faces,1);

    f_ids = find(arrayfun(@(x) ~isempty(x.BC_flux), face_struct));
    f_vals = [face_struct(f_ids).BC_flux];

    % for idx = 1:length(f_ids)
    %     f = f_ids(idx);
        A(f_ids, :) = 0;
        A(:, f_ids) = 0;
        A(f_ids, f_ids) = speye(length(f_ids));

        rhs_BC(f_ids) = f_vals;
    % end
end

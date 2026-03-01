function d_K = computeLocalDirichlet3D(c, cell_struct, face_struct, a, b, c3, d)

    face_ids = cell_struct(c).faces;
    signs = cell_struct(c).faces_orientation;
    nlf = length(face_ids);

    d_K = zeros(nlf,1);

    for j = 1:nlf
        f  = face_ids(j);
        xf = face_struct(f).center(1);
        yf = face_struct(f).center(2);
        zf = face_struct(f).center(3);

        p_face = a*xf + b*yf + c3*zf + d;
        d_K(j) = signs(j) * p_face;
    end
end
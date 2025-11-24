function [A, rhs_BC] = neumannBoundary_1d(A, cell_struct)
% Here we assume vertical fracture, so just need to modify two things:
% (1) Bottom vertex of first fracture segment
% (2) Top vertex of last fracture segment

    n_faces = length(cell_struct) + 1;
    rhs_BC = zeros(n_faces,1);

    A(:,1) = 0;
    A(1,:) = 0;
    A(1,1) = 1;

    A(:,n_faces) = 0;
    A(n_faces,:) = 0;
    A(n_faces,n_faces) = 1;

    rhs_BC(1) = 0;
    rhs_BC(n_faces) = 0;
   
end
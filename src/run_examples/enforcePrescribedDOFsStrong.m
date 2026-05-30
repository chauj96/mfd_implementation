function [A, b] = enforcePrescribedDOFsStrong( prescribedIdx, prescribedVal, A, b)
%enforcePrescribedDOFsStrong Strongly enforce prescribed DOFs values in A*x = b.
%   Enforces x(prescribedIdx) = prescribedVal by:
%     1) removing couplings to/from prescribed indices (zero rows and columns)
%     2) setting A(ii,ii) = 1 for prescribed indices
%     3) updating b so the reduced equations remain consistent

    nUnknowns = size(A, 1);

    % Optional checks
    assert(size(A,2) == nUnknowns, 'A must be square.');
    assert(numel(b) == nUnknowns, 'b must have length nUnknowns.');

    % Make prescribedVal a vector aligned with prescribedIdx
    if isscalar(prescribedVal)
        prescribedVal = repmat(prescribedVal, numel(prescribedIdx), 1);
    else
        prescribedVal = prescribedVal(:);
    end

    % Indicate which unknowns remain free
    isFree = true(nUnknowns, 1);
    isFree(prescribedIdx) = false;

    selectFree = spdiags(double(isFree), 0, nUnknowns, nUnknowns);
    selectPrescribed = spdiags(double(~isFree), 0, nUnknowns, nUnknowns);

    % Full-length vector with prescribed entries filled in
    xPrescribed = zeros(nUnknowns, 1);
    xPrescribed(prescribedIdx) = prescribedVal;

    % Update RHS to account for removing prescribed columns
    A_freeRows = selectFree * A;
    b = b - A_freeRows * xPrescribed;

    % Force equations at prescribed indices: x(i) = prescribedVal
    b(prescribedIdx) = prescribedVal;

    % Modify matrix: keep only free-free block plus identity on prescribed indices
    A = A_freeRows * selectFree + selectPrescribed;
end

% asec - arc secant function
function result = asec(this)
    % This function does not preserve the matrix sparsity: asec(0) ~= 0
    % But matlab's implementation returns a sparse result. So we provide
    % both implementations:

    % The cosine is a full matrix
    result = asec(full(this));

    % for matlab, cos of a sparse matrix is a sparse matrix
    if gemSparseLikeMatlab == 1
        result = sparse(result);
    end
end

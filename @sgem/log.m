% log - computes the natural logarithm
function result = log(this)
    % This function does not preserve the matrix sparsity: log(0) ~= 0
    % But matlab's implementation returns a sparse result. So we provide
    % both implementations:

    % The logarithm is a full matrix
    result = log(full(this));

    % for matlab, log of a sparse matrix is a sparse matrix
    if gemSparseLikeMatlab == 1
        result = sparse(result);
    end
end

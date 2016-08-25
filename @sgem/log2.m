% log2 - computes the logarithm in bases 2
function result = log2(this)
    % This function does not preserve the matrix sparsity: log2(0) ~= 0
    % But matlab's implementation returns a sparse result. So we provide
    % both implementations:

    % The logarithm is a full matrix
    result = log2(full(this));

    % for matlab, log of a sparse matrix is a sparse matrix
    if gemSparseLikeMatlab == 1
        result = sparse(result);
    end
end

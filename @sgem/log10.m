% log10 - computes the logarithm in bases 10
function result = log10(this)
    % This function does not preserve the matrix sparsity: log10(0) ~= 0
    % But matlab's implementation returns a sparse result. So we provide
    % both implementations:

    % The logarithm is a full matrix
    result = log10(full(this));

    % for matlab, log of a sparse matrix is a sparse matrix
    if gemSparseLikeMatlab == 1
        result = sparse(result);
    end
end

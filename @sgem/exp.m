% exp - computes the exponential function
function result = exp(this)
    % This function does not preserve the matrix sparsity: exp(0) ~= 0
    % But matlab's implementation returns a sparse result. So we provide
    % both implementations:

    % The exponential is a full matrix
    result = exp(full(this));

    % for matlab, exp of a sparse matrix is a sparse matrix
    if gemSparseLikeMatlab == 1
        result = sparse(result);
    end
end

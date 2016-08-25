% cot - cotangent function
function result = cot(this)
    % This function does not preserve the matrix sparsity: cot(0) ~= 0
    % But matlab's implementation returns a sparse result. So we provide
    % both implementations:

    % The cosine is a full matrix
    result = cot(full(this));

    % for matlab, cos of a sparse matrix is a sparse matrix
    if gemSparseLikeMatlab == 1
        result = sparse(result);
    end
end

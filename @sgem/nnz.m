% nnz - number of non-zero elements in the matrix
function result = nnz(this)
    % We call the nnz procedure.
    result = sgem_mex('nnz', this.objectIdentifier);
end

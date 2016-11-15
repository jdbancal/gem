% nnz - number of non-zero elements in the matrix
function result = nnz(this)
    % We obtain the result from the find method
    result = length(find(this));
end

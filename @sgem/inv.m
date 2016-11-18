% inv - matrix inverse (computed by LU decomposition)
function result = inv(this)
    % Inverse matrices are only defined for square matrices
    if size(this,1) ~= size(this,2)
        error('Matrix must be square in sgem::inv');
    end
    
    % We first check whether the rank of the matrix is maximal
    if rank(this) < size(this,1)
        warning('Matrix is singular to working precision');
        result = Inf*gem(ones(size(this)));
        return;
    end
    
    % We call the inv procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = sgem_mex('inv', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

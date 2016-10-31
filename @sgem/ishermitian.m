% ishermitian - tells whether the matrix is hermitian, i.e. x == x'
function result = ishermitian(this)
    % We call the ishermitian procedure
    result = logical(sgem_mex('ishermitian', this.objectIdentifier));
end

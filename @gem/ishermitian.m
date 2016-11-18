% ishermitian - tells whether the matrix is hermitian, i.e. x == x'
function result = ishermitian(this)
    % We call the ishermitian procedure
    result = logical(gem_mex('ishermitian', this.objectIdentifier));
end

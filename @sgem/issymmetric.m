% issymmetric - tells whether the matrix is symmetric, i.e. x == x.'
function result = issymmetric(this)
    % We call the ishermitian procedure
    result = logical(sgem_mex('issymmetric', this.objectIdentifier));
end

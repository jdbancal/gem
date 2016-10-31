% issymmetric - tells whether the matrix is symmetric, i.e. x == x.'
function result = issymmetric(this)
    % We call the issymmetric procedure
    result = logical(gem_mex('issymmetric', this.objectIdentifier));
end

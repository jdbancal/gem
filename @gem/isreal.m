% isreal - tells whether the matrix is real
function result = isreal(this)
    % We call the isreal procedure
    result = gem_mex('isreal', this.objectIdentifier);
end

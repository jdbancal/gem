% isreal - tells whether the matrix is real
function result = isreal(this)
    % We call the isreal procedure
    result = sgem_mex('isreal', this.objectIdentifier);
end

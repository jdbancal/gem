% rank - number of linearly independent lines in the matrix
function result = rank(this)
    % We call the rank procedure.
    result = sgem_mex('rank', this.objectIdentifier);
end

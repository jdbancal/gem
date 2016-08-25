% sparse - converts the matrix into sparse format
function result = sparse(this)
    % We call the sgem constructor with the current object
    result = sgem(this);
end

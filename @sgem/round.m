% round - rounds towards the nearest integer
function result = round(this)
    % We call the rounding procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = sgem_mex('round', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

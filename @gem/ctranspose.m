% ctranspose - computes the hermitian conjugate
function result = ctranspose(this)
    % We call the ctranspose procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('ctranspose', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

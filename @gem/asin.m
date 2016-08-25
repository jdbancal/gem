% asin - arc sine function
function result = asin(this)
    % We call the asin procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('asin', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end
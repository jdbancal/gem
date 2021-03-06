% tan - tangent function
function result = tan(this)
    % We call the tan procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('tan', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

% inv - matrix inverse (computed by LU decomposition)
function result = inv(this)
    % Inverse matrices are only defined for square matrices
    if size(this,1) ~= size(this,2)
        error('Matrix must be square in gem::inv');
    end

    % We call the inv procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('inv', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

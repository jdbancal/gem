% sqrt - element-wise square root
function result = sqrt(this)
    % A non-optimized but simple way to compute the square root
%result = this.^(1/gem(2));

    % We call the sqrt procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('sqrt', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

% cbrt - element-wise cubic root (for positive real numbers)
function result = cbrt(this)
    % First we check that the matrix only contains positive real numbers
    if (~isreal(this)) || (min(min(this)) < 0)
        % Then we cannot use the quick computation. We use a generic one instead
        result = this.^(1/gem(3));
        return;
    end

    % We call the cbrt procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('cbrt', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

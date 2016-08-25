% cbrt - element-wise cubic root (for positive real numbers)
function result = cbrt(this)
    % First we check that the matrix only contains positive real numbers
    if (~isreal(this)) || (min(min(full(this))) < 0) % REMOVE THE FULL HERE ONCE THE MIN FUNCTION IS IMPLEMENTED FOR SGEM OBJECTS
        % Then we cannot use the quick computation. We use a generic one instead
        result = this.^(1/gem(3));
        return;
    end

    % We call the cbrt procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = sgem_mex('cbrt', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

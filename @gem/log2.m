% log2 - computes the logarithm in bases 2
function result = log2(this)
    % We call the log procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('log', this.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
    
    % We put the result in the correct basis
    % Warning : Here we use the default current precision to construct the 
    %   rescaling coefficient... This precision may be smaller than that of 
    %   the numbers in the current matrix, hence resulting in an result with 
    %   bad precision...
    result = result/log(gem('2'));
end

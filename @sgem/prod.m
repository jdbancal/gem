% prod - Product of elements
%
% supported formats :
%   prod(a) : column-wise product
%   prod(a, dim) : product along the dimension dim
function result = prod(this, dim)
    % This function can involve up to two arguments
    if (nargin >= 2) && (iscell(dim) || (numel(dim) > 1) || (dim ~= 1 && dim ~= 2))
        error('Unexpected arguments in sgem::prod');
    end

    if nargin < 2
        if size(this,1) ~= 1
            dim = 1;
        else
            dim = 2;
        end
    end
    
    % Now we call the element-wise minimum procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    switch dim
        case 1
            newObjectIdentifier = sgem_mex('colProd', this.objectIdentifier);
        case 2
            newObjectIdentifier = sgem_mex('rowProd', this.objectIdentifier);
        otherwise
            error('Unexpected argument in sgem::prod');
    end
    
    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

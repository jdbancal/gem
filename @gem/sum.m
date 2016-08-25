% sum - Sum of elements
%
% supported formats :
%   sum(a) : column-wise sum
%   sum(a, dim) : sum along the dimension dim
function result = sum(this, dim)
    % This function can involve up to two arguments
    if (nargin >= 2) && (iscell(dim) || (numel(dim) > 1) || (dim ~= 1 && dim ~= 2))
        error('Unexpected arguments in gem::sum');
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
            newObjectIdentifier = gem_mex('colSum', this.objectIdentifier);
        case 2
            newObjectIdentifier = gem_mex('rowSum', this.objectIdentifier);
        otherwise
            error('Unexpected argument in gem::sum');
    end
    
    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

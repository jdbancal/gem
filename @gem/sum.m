% sum - Sum of elements
%
% supported formats :
%   sum(a) : column-wise sum
%   sum(a, dim) : sum along the dimension dim
%   sum(a, dim, type) : sum along the dimension dim and type can be any of
%       the following:
%        - double : the output will be recast in double type
%        - native or default : the output is of gem type
function result = sum(this, dim, type)
    % This function can involve up to three arguments
    if (nargin < 3)
        type = 'native';
    else
        if ~isequal(type, 'double') && ~isequal(type, 'native') && ~isequal(type, 'default')
            error('Unsupported argument in gem::sum');
        end
    end
    
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
    
    % If the user requested double output, we transform the output
    % accordingly
    if isequal(type, 'double')
        result = double(result);
    end
end

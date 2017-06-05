% sum - Sum of elements
%
% supported formats :
%   sum(a) : column-wise sum
%   sum(a, dim) : sum along the dimension dim
%   sum(a, dim, type) : sum along the dimension dim and type can be any of
%       the following:
%        - double : the output will be recast in double type
%        - native or default : the output is of sgem type
function result = sum(this, dim, type)
    % This function can involve up to three arguments
    if (nargin < 3)
        type = 'native';
    else
        if ~isequal(type, 'double') && ~isequal(type, 'native') && ~isequal(type, 'default')
            error('Unsupported argument in sgem::sum');
        end
    end
    
    if (nargin >= 2) && (iscell(dim) || (numel(dim) > 1) || (dim ~= 1 && dim ~= 2))
        error('Unexpected arguments in sgem::sum');
    end

    if nargin < 2
        if size(this,1) ~= 1
            dim = 1;
        else
            dim = 2;
        end
    end
    
    % Now we compute the sum by multiplying with a unit vector.
    switch dim
        case 1
            result = sparse(ones(1,size(this,1)))*this;
        case 2
            result = this*sparse(ones(size(this,2),1));
        otherwise
            error('Unexpected argument in sgem::sum');
    end
    
    % If the user requested double output, we transform the output
    % accordingly
    if isequal(type, 'double')
        result = double(result);
    end
end

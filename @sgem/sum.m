% sum - Sum of elements
%
% supported formats :
%   sum(a) : column-wise sum
%   sum(a, dim) : sum along the dimension dim
function result = sum(this, dim)
    % This function can involve up to two arguments
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
    
end

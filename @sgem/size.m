% size of the matrix
function [result1 result2] = size(this, dim)
    if (nargin >= 2) && (dim ~= 1) && (dim ~= 2)
        error('Wrong dimension in sgem::size');
    end

    result1 = sgem_mex('size', this.objectIdentifier);
    
    if nargin >= 2
        result1 = result1(dim);
    elseif nargout == 2
        result2 = result1(2);
        result1 = result1(1);
    end
end

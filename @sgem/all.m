% all - True if all elements are nonzero
%
% supported formats:
%   all(x)
%   all(x,dim)
function result = all(this, dim)

if (nargin < 2) && (min(size(this)) > 1)
    % For matrices, we act by default on the columns
    dim = 1;
elseif nargin < 2
    result = (sum(this~=0) == size(this,2));
    return;
end

% We look along a given dimension
if numel(dim) ~= 1 || (round(dim) ~= dim) || (dim < 1)
    error('Dim must be a positive integer');
end
result = (sum(this~=0,dim) == size(this,dim));

end

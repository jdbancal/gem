% any - True for numbers, ignores zeros and NaN elements
%
% supported formats:
%   any(x)
%   any(x,dim)
function result = any(this, dim)

if (nargin < 2) && (min(size(this)) > 1)
    % For matrices, we act by default on the columns
    dim = 1;
elseif nargin < 2
    result = logical(sum((this ~= 0) - isnan(this)));
    return;
end

% We look along a given dimension
if numel(dim) ~= 1 || (round(dim) ~= dim) || (dim < 1)
    error('Dim must be a positive integer');
end
result = logical(sum((this ~= 0) - isnan(this),dim));

end

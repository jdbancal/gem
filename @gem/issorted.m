% issorted - checks if a matrix is sorted
%
% supported formats:
%   issorted(vector)
%   issorted(matrix, 'rows')
function result = issorted(this, option)
    if nargin < 2
        if min(size(this)) > 1
            error('Input must be a vector, or option ''rows'' ust be used');
        end
        result = (this == sort(this));
    else
        if ~isequal(option, 'rows')
            error('Unknown option');
        end
        result = (this == sortrows(this));
    end
end

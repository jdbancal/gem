% diff - Sum of elements
%
% supported formats :
%   diff(a) : column-wise difference between successive elements
%   diff(a, n) : n-th difference
%   diff(a, n, dim) : n-th difference along dimension dim
function result = diff(this, n, dim)
    % Input management
    if nargin < 3
        if size(this,1) ~= 1
            dim = 1;
        else
            dim = 2;
        end
    else
        if (numel(dim) ~= 1) || (dim < 1) || (round(dim) ~= dim)
            error('Dim must be a positive integer');
        end
    end
    
    if nargin < 2
        n = 1;
    else
        if (numel(n) ~= 1) || (n < 1) || (round(n) ~= n)
            error('N must be an integer');
            % Here we could support more than matlab by allowing 
            % differences of order 0 (and returns itself)
        end
    end
    
    % Producing the result
    if numel(this) == 1
        result = gem([]);
        return;
    end
    
    if n == 0
        result = this;
        return;
    end
    
    if n > 1
        this = diff(this, n-1, dim);
    end
    
    if dim == 1
        subs1.type = '()';
        subs1.subs = {[2:size(this,1)] ':'};
        subs2.type = '()';
        subs2.subs = {[1:size(this,1)-1] ':'};
        result = subsref(this, subs1) - subsref(this, subs2);
    else
        subs1.type = '()';
        subs1.subs = {':' [2:size(this,2)]};
        subs2.type = '()';
        subs2.subs = {':' [1:size(this,2)-1]};
        result = subsref(this, subs1) - subsref(this, subs2);
    end
end

% sort - sorts components
%
% supported formats :
%   sort(a) : sorts out elements of a along the first dimension in
%     ascending order
%   sort(a, dim) : sorts out elements of a along the first dimension dim
%   sort(a, mode) : mode is either 'ascend' or 'descend'
%   sort(a, dim, mode) : specifies both the dimension and the mode
%   [b I] = sort(a) : returns the indices I such that b = a(I)
function [result I] = sort(this, dim, mode)
    % Input management
    if nargin < 3
        mode = 'ascend';
    end
    
    if nargin < 2
        if size(this,1) ~= 1
            dim = 1;
        else
            dim = 2;
        end
    elseif ~isequal(class(dim), 'double')
        dim = double(dim);
    end
    
    if ischar(dim)
        mode = dim;
        dim = 1;
    end
    
    if ~isequal(mode, 'ascend') && ~isequal(mode, 'descend')
        error('Sorting direction must be either ''ascend'' of ''descend''');
    end
    
    if (numel(dim) ~= 1) || (dim < 1) || (round(dim) ~= dim)
        error('Dim must be a positive integer');
    end

    % If there is nothing to sort
    if dim >= 3
        result = this;
        if nargout > 1
            if size(this,1) ~= 1
                I = [1:size(this,1)]'*ones(1,size(this,2));
            else
                I = [1:size(this,2)];
            end
        end
        return;
    end
    
    % Now we call the appropriate sorting method
    [newObjectIdentifier I] = gem_mex('sort', this.objectIdentifier, dim-1, double(isequal(mode, 'descend')));

    % Indices in matlab start from 1
    I = I+1;
    
    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

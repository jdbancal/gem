% sort - sorts components
%
% supported formats :
%   sort(a) : sorts out elements of a along the first dimension in
%     ascending order
%   sort(a, dim) : sorts out elements of a along the first dimension dim
%   sort(a, mode) : mode is either 'ascend' or 'descend'
%   sort(a, dim, mode) : specifies both the dimension and the mode
%   [b I] = sort(a) : returns the indices I such that b = a(I)
%     (warning : indices are a full matrix)
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
    [newObjectIdentifier Icell nbNegatives] = sgem_mex('sort', this.objectIdentifier, dim-1, double(isequal(mode, 'descend')));

    % Indices are in a compressed form; we process them if necessary
    if nargout > 1
        I = zeros(size(this));
        if dim == 2
            I = I';
        end
        
        for k = 1:size(I,2)
            Ik = Icell{k}+1;
            
            if isequal(mode, 'ascend')
                startBlock = nbNegatives(k);
                endBlock = length(Ik) - nbNegatives(k);
            else
                startBlock = length(Ik) - nbNegatives(k);
                endBlock = nbNegatives(k);
            end
            
            I(1:startBlock, k) = Ik(1:startBlock);
            I(startBlock+1:end-endBlock, k) = setdiff(1:size(I,1),Ik);
            I(end-endBlock+1:end, k) = Ik(startBlock+1:end);
        end

        if dim == 2
            I = I';
        end
    end
    
    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

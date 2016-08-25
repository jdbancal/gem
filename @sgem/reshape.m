% reshape - reshapes the matrix (keeping the number of elements constant)
function result = reshape(this, varargin)

    if length(varargin) == 1
        if numel(varargin{1}) ~= 2
            error('Wrong number dimensions');
        end
        newDim = varargin{1};
    elseif length(varargin) == 2
        if isempty(varargin{1}) && isempty(varargin{2})
            error('Not enough dimensions specified');
        end
        if ~isempty(varargin{1})
            newDim(1) = varargin{1};
        else
            newDim(1) = numel(this)/varargin{2};
        end
        if ~isempty(varargin{2})
            newDim(2) = varargin{2};
        else
            newDim(2) = numel(this)/newDim(1);
        end
    else
        error('Inappropriate call to reshape');
    end

    if prod(newDim) ~= numel(this)
        error('Reshape cannot change the number of elements');
    end
    
    % We call the abs procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = sgem_mex('reshape', this.objectIdentifier, newDim(1), newDim(2));

    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

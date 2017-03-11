% rdivide - element-wise right division A./B
function result = rdivide(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in sgem::rdivide');
    end

    % If the second type is more elaborated than just a number, we let the
    % corresponding class take care of performing the operation
    if ~isnumeric(varargin{1})
        result = ldivide(varargin{1}, this);
        return;
    end
    
    % This function preserves the sparsity only if A is sparse and B is a scalar.
    if issparse(this) && (numel(varargin{1}) == 1)
        % We check that both objects are of type sgem,
        % otherwise we convert them
        if ~isequal(class(this), 'sgem')
            this = sgem(this);
        elseif ~isequal(class(varargin{1}), 'sgem')
            varargin{1} = sgem(varargin{1});
        end

        % We call the rdivide procedure. Since the function creates a
        % new object with the result, we keep the corresponding handle...
        newObjectIdentifier = sgem_mex('rdivide', this.objectIdentifier, varargin{1}.objectIdentifier);

        % ...  and create a new matlab object to keep this handle
        result = sgem('encapsulate', newObjectIdentifier);
    else
        % If the result is not sparse, we compute this function on the full variables
        result = rdivide(full(this), full(varargin{1}));

        % For matlab, division of a sparse matrix is always a sparse matrix
        if gemSparseLikeMatlab == 1
            result = sparse(result);
        end
    end
end

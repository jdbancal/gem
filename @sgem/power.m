% power - element-wise power A.^B
function result = power(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in sgem::power');
    end

    % We need to check that the operation is possible (the c++
    % library might give bad errors otherwise). So we request the
    % dimensions of each matrix
    size1 = size(this);
    size2 = size(varargin{1});

    if (~isequal(size1, size2)) && (prod(size1) ~= 1) && (prod(size2) ~= 1)
        error('Incompatible size for element-wise matrix power');
    end

    % The result of a sparse power is sparse if the exponent is > 0
    % (note that matlab uses a slightly different convention)
    if numel(this) >= numel(varargin{1}) && isreal(varargin{1}) && (min(min(full(varargin{1}))) > 0) % REMOVE THE FULL ONCE MIN AND > ARE IMPLEMENTED FOR SPARSE MATRICES
        % To get a sparse result, the current matrix must be sparse
        if ~isequal(class(varargin{1}), 'sgem')
            varargin{1} = sgem(varargin{1});
        end
        % The exponent must be in full form (it contains no zero)
        if ~isequal(class(varargin{1}), 'gem')
            varargin{1} = gem(varargin{1});
        end

        % Now we call the power procedure. Since the function creates a
        % new object with the result, we keep the corresponding handle...
        newObjectIdentifier = sgem_mex('power', this.objectIdentifier, varargin{1}.objectIdentifier);

        % ...  and create a new matlab object to keep this handle
        result = sgem('encapsulate', newObjectIdentifier);
    else
        % In this case the result is a full matrix, so we translate
        % the matrices to full form (and gem)
        if ~isequal(class(this), 'gem')
            this = gem(this);
        elseif ~isequal(class(varargin{1}), 'gem')
            varargin{1} = gem(varargin{1});
        end
        % ... and call the function for dense matrices
        result = this.^varargin{1};

        % For matlab, power of a sparse matrix is always sparse
        if gemSparseLikeMatlab == 1
            result = sparse(result);
        end
    end

end

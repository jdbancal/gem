% times - element-wise product A.*B
function result = times(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in gem::times');
    end

    % If the second type is more elaborated than just a number, we let the
    % corresponding class take care of performing the multiplication
    if ~isnumeric(varargin{1})
        result = times(varargin{1}, this);
        return;
    end
    
    % We need to check that the operation is possible (the c++
    % library might give bad errors otherwise). So we request the
    % dimensions of each matrix
    size1 = size(this);
    size2 = size(varargin{1});

    if (~isequal(size1, size2)) && (prod(size1) ~= 1) && (prod(size2) ~= 1)
        error('Incompatible size for element-wise matrix multiplication');
    end

    % Now we also check the type of both objects
    % First we make sure that both object are either full of sparse gems
    if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
        this = gemify(this);
    end
    if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = gemify(varargin{1});
    end

    % We only implement multiplication by a scalar when both the matrix and
    % the scalar are full, or when they are both sparse. Element-wise
    % multiplication between two 'sparse' matrices anyway produces a sparse
    % result.
    if (prod(size1) == 1) && ~isequal(class(this), class(varargin{1}))
        if isequal(class(this), 'gem')
            this = sparse(this);
        else
            this = full(this);
        end
    end
    if (prod(size2) == 1) && ~isequal(class(this), class(varargin{1}))
        if isequal(class(varargin{1}), 'gem')
            varargin{1} = sparse(varargin{1});
        else
            varargin{1} = full(varargin{1});
        end
    end
    
    % Now we call the appropriate multiplication procedure, and store the
    % result in the adequate object.
    if isequal(class(this), 'gem')
        if isequal(class(varargin{1}), 'gem')
            newObjectIdentifier = gem_mex('times', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = gem('encapsulate', newObjectIdentifier);
        else
            newObjectIdentifier = gem_mex('times_fs', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = sgem('encapsulate', newObjectIdentifier);
        end
    else
        % A priori we should not arrive here... but just in case
        if isequal(class(varargin{1}), 'gem')
            newObjectIdentifier = sgem_mex('times_sf', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = sgem('encapsulate', newObjectIdentifier);
        else
            newObjectIdentifier = sgem_mex('times', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = sgem('encapsulate', newObjectIdentifier);
        end
    end
    
end


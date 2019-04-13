% mtimes - multiplies two matrices
function result = mtimes(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in gem::plus');
    end

    % If the second type is more elaborated than just a number, we let the
    % corresponding class take care of performing the multiplication
    if ~isnumeric(varargin{1})
        result = mtimes(varargin{1}.', this.').';
        return;
    end
    
    % We need to check that the operation is possible (the c++
    % library might give bad errors otherwise). So we request the
    % dimensions of each matrix
    size1 = size(this);
    size2 = size(varargin{1});

    % If one number is a scalar, we reinterpret the multiplication as
    % an element-wise operation
    if (prod(size1) == 1) || (prod(size2) == 1)
        result = this.*varargin{1};
        return;
    end
    
    % We check that the sizes are compatible
    if size1(2) ~= size2(1)
        error('Incompatible size for matrix multiplication');
    end

    % Now we also check the type of both objects
    % First we make sure that both object are either full of sparse gems
    if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
        this = gemify(this);
    end
    if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = gemify(varargin{1});
    end
    
    % Now we call the appropriate multiplication procedure, and store the
    % result in the adequate object.
    if isequal(class(this), 'gem')
        if isequal(class(varargin{1}), 'gem')
            newObjectIdentifier = gem_mex('mtimes', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = gem('encapsulate', newObjectIdentifier);
        else
            newObjectIdentifier = gem_mex('mtimes_fs', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = gem('encapsulate', newObjectIdentifier);
        end
    else
        % A priori we should not arrive here... but just in case
        if isequal(class(varargin{1}), 'gem')
            newObjectIdentifier = sgem_mex('mtimes_sf', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = gem('encapsulate', newObjectIdentifier);
        else
            newObjectIdentifier = sgem_mex('mtimes', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = sgem('encapsulate', newObjectIdentifier);
        end
    end
    
end

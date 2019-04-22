% mldivide
%
% Only simplest cases are supported at the moment
function result = mldivide(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in sgem::mldivide');
    end
    
    % If the second type is more elaborated than just a number, we let the
    % corresponding class take care of performing the operation
    if ~isnumeric(varargin{1})
        result = mrdivide(varargin{1}', this')';
        return;
    end
    
    % We check that the denominator is a scalar
    if numel(this) == 1
        result = varargin{1}./this;
        return;
    end
    
    % We check whether the dimensions match
    if size(this,1) ~= size(varargin{1},1)
        error('Size mismatch in sgem::mldivide');
    end
    
    % We check if the denominator is singular
    if rank(this) < size(this,1)
        error('Matrix is singular in sgem::mldivide');
    end
    
    % We check the conditioning number of the matrix
    rcond = cond(this);
    if rcond > 10^gemWorkingPrecision
        warning(['Matrix is close to singular. Result may be inaccurate. RCOND = ', toStrings(rcond)]);
    end
    
    % Now we also check the type of both objects
    % First we make sure that both object are either full of sparse gems
    if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
        this = gemify(this);
    end
    if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = gemify(varargin{1});
    end
    % We only support solving a dense system with a dense solution
    if isequal(class(this), 'gem') && isequal(class(varargin{1}), 'sgem')
        varargin{1} = full(varargin{1});
    end
    
    % Now we call the appropriate division procedure, and store the
    % result in the adequate object.
    if isequal(class(this), 'gem')
        newObjectIdentifier = gem_mex('mldivide', this.objectIdentifier, varargin{1}.objectIdentifier);
        result = gem('encapsulate', newObjectIdentifier);
    else
        % A priori we should not arrive here... but just in case
        if isequal(class(varargin{1}), 'gem')
            newObjectIdentifier = sgem_mex('mldivide_sf', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = gem('encapsulate', newObjectIdentifier);
        else
            newObjectIdentifier = sgem_mex('mldivide', this.objectIdentifier, varargin{1}.objectIdentifier);
            result = sgem('encapsulate', newObjectIdentifier);
        end
    end

end

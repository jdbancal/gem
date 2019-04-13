% lt - lesser than
% 
% For complex numbers, this function compares only the real parts
function result = lt(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in sgem::lt');
    end

    % We need to check that the operation is possible (the c++
    % library might give bad errors otherwise). So we request the
    % dimensions of each matrix
    size1 = size(this);
    size2 = size(varargin{1});

    if (~isequal(size1, size2)) && (prod(size1) ~= 1) && (prod(size2) ~= 1)
        error('Incompatible size for lt comparison');
    end

    % We make sure both objects are gem objects
    if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
        this = gemify(this);
    end
    if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = gemify(varargin{1});
    end

    % comparison of a full matrix with a sparse one gives a sparse result
    if issparse(this) && (~issparse(varargin{1}))
        varargin{1} = sparse(varargin{1});
    end
    if ~issparse(this) && (issparse(varargin{1}))
        this = sparse(this);
    end

    % If one object is sparse and the comparison involves a scalar, we
    % check whether the comparison leaves zero elements in the matrix 
    % invariant. If this is not the case, the output will be full, so we
    % convert the objects to full
    convertBackForMatlabStyle = 0;
    if (issparse(this) && (prod(size2) == 1) && (0 < full(varargin{1}))) ...
            || (issparse(varargin{1}) && (prod(size1) == 1) && (full(this) < 0))
        this = full(this);
        varargin{1} = full(varargin{1});
        convertBackForMatlabStyle = 1;
    end
    
    % Now we can compute the comparison matrix
    if isequal(class(this), 'gem') && isequal(class(varargin{1}), 'gem')
        result = logical(gem_mex('lt', this.objectIdentifier, varargin{1}.objectIdentifier));
    else
        result = logical(sgem_mex('lt', this.objectIdentifier, varargin{1}.objectIdentifier));
    end
    
    % For matlab, lt when one object is sparse produces a sparse result
    if (gemSparseLikeMatlab == 1) && (convertBackForMatlabStyle == 1)
        result = sparse(result);
    end
end

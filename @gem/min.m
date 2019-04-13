% min - smallest component
%
% supported formats :
%   [Y I] = min(a) : column-wise minimum (or line-wise for a line vector)
%   [Y I] = min(a,[],1) : column-wise minimum
%   [Y I] = min(a,[],2) : line-wise minimum
%   Y = min(a,b) : element-wise minimum (either a or b can also be scalars)
function [Y I] = min(this, varargin)
    % This function can involve up to three arguments
    if length(varargin) > 2
        error('Wrong number of arguments in gem::min');
    end

    % We check in which way the function is called
    if (length(varargin) == 2) && (~isempty(varargin{1}) || (~isequal(varargin{2},1) && ~isequal(varargin{2},2)))
        error('Wrong arguments for gem::min');
    end
    
    if (length(varargin) == 1) 
        % We need to check that the operation is possible (the c++
        % library might give bad errors otherwise). So we request the
        % dimensions of each matrix
        size1 = size(this);
        size2 = size(varargin{1});

        if (~isequal(size1, size2)) && (prod(size1) ~= 1) && (prod(size2) ~= 1)
            error('Incompatible sizes for element-wise minimum');
        end
        
        % We also check that no second argument is expected
        if nargout > 1
            error('Element-wise minimum returns only one outcome');
        end
    end
    
    
    %% If we reach here, the arguments must be good    
    if length(varargin) ~= 1
        size1 = size(this);
        % Now we call the column or line-wise minimum procedure. Since the function creates a
        % new object with the result, we keep the corresponding handle...
        if (isempty(varargin) && (size1(1) > 1)) || ((length(varargin) == 2) && (varargin{2} == 1))
            [newObjectIdentifier I] = gem_mex('colMin', this.objectIdentifier);
            I = I'+1;
        else
            [newObjectIdentifier I] = gem_mex('rowMin', this.objectIdentifier);
            I = I+1;
        end
        % ...  and create a new matlab object to keep this handle
        Y = gem('encapsulate', newObjectIdentifier);
    else
        % Here we compute the minimum between two objects
        
        % For later, we remember if one of the input is sparse
        oneInputSparse = max(1, issparse(this) + issparse(varargin{1}));

        % Let's we make sure that both object are either full of sparse gems
        if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
            this = gemify(this);
        end
        if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
            varargin{1} = gemify(varargin{1});
        end
        
        size1 = size(this);
        size2 = size(varargin{1});
        
        % We only implement minimum with a scalar when both the matrix and
        % the scalar are full, or when both are sparse. Minimum between two
        % matrices is sparse by default.
        if (prod(size1) == 1) && (prod(size1) == 1) && ~isequal(class(this), class(varargin{1}))
            this = sparse(this);
            varargin{1} = sparse(varargin{1});
        elseif (prod(size1) == 1) && ~isequal(class(this), class(varargin{1}))
            if isequal(class(this), 'gem')
                this = sparse(this);
            else
                this = full(this);
            end
        elseif (prod(size2) == 1) && ~isequal(class(this), class(varargin{1}))
            if isequal(class(varargin{1}), 'gem')
                varargin{1} = sparse(varargin{1});
            else
                varargin{1} = full(varargin{1});
            end
        elseif ~isequal(class(this), class(varargin{1}))
            % Without further garantee, minimum between sparse and full
            % matrices could be full, but it may not be, so we produce a 
            % sparse output
            this = sparse(this);
            varargin{1} = sparse(varargin{1});
        end

        % Real minimum with a negative scalar is always full, so we 
        % produce a full output in this case
        if isreal(this) && isreal(varargin{1}) && (((numel(this) == 1) && (numel(varargin{1}) ~= 1) && (this < 0)) || ((numel(this) ~= 1) && (numel(varargin{1}) == 1) && (varargin{1} < 0)))
            % In this case, we compute the minimum of the full matrices
            this = full(this);
            varargin{1} = full(varargin{1});
        end
        
        if ~issparse(this)
            % Now we call the element-wise minimum procedure. Since the function creates a
            % new object with the result, we keep the corresponding handle...
            newObjectIdentifier = gem_mex('ewMin', this.objectIdentifier, varargin{1}.objectIdentifier);
            % ...  and create a new matlab object to keep this handle
            Y = gem('encapsulate', newObjectIdentifier);
        else
            % Now we call the element-wise minimum procedure. Since the function creates a
            % new object with the result, we keep the corresponding handle...
            newObjectIdentifier = sgem_mex('ewMin', this.objectIdentifier, varargin{1}.objectIdentifier);
            % ...  and create a new matlab object to keep this handle
            Y = sgem('encapsulate', newObjectIdentifier);
        end

        % For matlab, min when one object is sparse produces a sparse result
        if (gemSparseLikeMatlab == 1) && (oneInputSparse == 1) && (issparse(Y) == 0)
            Y = sparse(Y);
        end
    end

end

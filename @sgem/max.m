% max - largest component
%
% supported formats :
%   [Y I] = max(a) : column-wise maximum (or line-wise for a line vector)
%   [Y I] = max(a,[],1) : column-wise maximum
%   [Y I] = max(a,[],2) : line-wise maximum
%   Y = max(a,b) : element-wise maximum (either a or b can also be scalars)
function [Y I] = max(this, varargin)
    % This function can involve up to three arguments
    if length(varargin) > 2
        error('Wrong number of arguments in sgem::max');
    end

    % We check in which way the function is called
    if (length(varargin) == 2) && (~isempty(varargin{1}) || (~isequal(varargin{2},1) && ~isequal(varargin{2},2)))
        error('Wrong arguments for sgem::max');
    end
    
    if (length(varargin) == 1) 
        % We need to check that the operation is possible (the c++
        % library might give bad errors otherwise). So we request the
        % dimensions of each matrix
        size1 = size(this);
        size2 = size(varargin{1});

        if (~isequal(size1, size2)) && (prod(size1) ~= 1) && (prod(size2) ~= 1)
            error('Incompatible sizes for element-wise maximum');
        end
        
        % We also check that no second argument is expected
        if nargout > 1
            error('Element-wise maximum returns only one outcome');
        end
    end
    
    
    %% If we reach here, the arguments must be good    
    if length(varargin) ~= 1
        size1 = size(this);
        % Now we call the column or line-wise maximum procedure. Since the function creates a
        % new object with the result, we keep the corresponding handle...
        if (isempty(varargin) && (size1(1) > 1)) || ((length(varargin) == 2) && (varargin{2} == 1))
            [newObjectIdentifier I] = sgem_mex('colMax', this.objectIdentifier);
            I = I'+1;
        else
            [newObjectIdentifier I] = sgem_mex('rowMax', this.objectIdentifier);
            I = I+1;
        end
        % ...  and create a new matlab object to keep this handle
        Y = sgem('encapsulate', newObjectIdentifier);
    else
        % Here we compute the maximum between two objects
        
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
        
        % We only implement maximum with a scalar when both the matrix and
        % the scalar are full, or when both are sparse. Maximum between two
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
            % Without further garantee, maximum between sparse and full
            % matrices could be full, but it may not be, so we produce a 
            % sparse output
            this = sparse(this);
            varargin{1} = sparse(varargin{1});
        end

        % Maximum with a positive scalar, as well as complex maximum with a
        % non-zero scalar, is always full, so we produce a full output in 
        % this case
        if (numel(this) == 1) && (numel(varargin{1}) ~= 1)
            if (~isreal(this) || ~isreal(varargin{1})) && (abs(this) == 0)
                % Then the max with zero is itself -> nothing to do, and we
                % keep sparse matrices sparse
                Y = varargin{1};
                return;
            elseif real(this) > 0
                % In this case, we compute the minimum of the full matrices
                this = full(this);
                varargin{1} = full(varargin{1});
            end
        end
        if (numel(this) ~= 1) && (numel(varargin{1}) == 1)
            if (~isreal(this) || ~isreal(varargin{1})) && (abs(varargin{1}) == 0)
                % Then the max with zero is itself -> nothing to do, and we
                % keep sparse matrices sparse
                Y = this;
                return;
            elseif real(varargin{1}) > 0
                % In this case, we compute the minimum of the full matrices
                this = full(this);
                varargin{1} = full(varargin{1});
            end
        end
        
        if ~issparse(this)
            % Now we call the element-wise minimum procedure. Since the function creates a
            % new object with the result, we keep the corresponding handle...
            newObjectIdentifier = gem_mex('ewMax', this.objectIdentifier, varargin{1}.objectIdentifier);
            % ...  and create a new matlab object to keep this handle
            Y = gem('encapsulate', newObjectIdentifier);
        else
            % Now we call the element-wise minimum procedure. Since the function creates a
            % new object with the result, we keep the corresponding handle...
            newObjectIdentifier = sgem_mex('ewMax', this.objectIdentifier, varargin{1}.objectIdentifier);
            % ...  and create a new matlab object to keep this handle
            Y = sgem('encapsulate', newObjectIdentifier);
        end

        % For matlab, max when one object is sparse produces a sparse result
        if (gemSparseLikeMatlab == 1) && (oneInputSparse == 1) && (issparse(Y) == 0)
            Y = sparse(Y);
        end
    end
    
end

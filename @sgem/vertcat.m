% vertcat - Vertical concatenation
function result = vertcat(this, varargin)
    % If the concatenation involves more than 2 objects, we deal with each pair
    % independently
    if length(varargin) > 1
        result = vertcat(this, varargin{1});
        for i = 2:length(varargin)
            result = vertcat(result, varargin{i});
        end
        return;
    end
    % Here there are only two elements to concatenate

    % If one element is empty, the concatenation is the other element
    if isempty(this)
        result = varargin{1};
        return;
    elseif isempty(varargin{1})
        result = this;
        return;
    end

    % We check that the two objects have a compatible size
    s1 = size(this);
    s2 = size(varargin{1});
    if s1(2) ~= s2(2)
        error('Incompatible sizes for sgem:vertcat');
    end

    % We make sure both objects are gem objects
    if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
        this = gemify(this);
    end
    if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = gemify(varargin{1});
    end

    % Concatenation of a full matrix with a sparse one gives a sparse result
    if issparse(this) && (~issparse(varargin{1}))
        varargin{1} = sparse(varargin{1});
    end
    if ~issparse(this) && (issparse(varargin{1}))
        this = sparse(this);
    end

    %% Now we can concatenate the two objects
    if isequal(class(this), 'gem') && isequal(class(varargin{1}), 'gem')
        newObjectIdentifier = gem_mex('vertcat', objectIdentifier(this), objectIdentifier(varargin{1}));
        result = gem('encapsulate', newObjectIdentifier);
    else
        newObjectIdentifier = sgem_mex('vertcat', objectIdentifier(this), objectIdentifier(varargin{1}));
        result = sgem('encapsulate', newObjectIdentifier);
    end
end

% vertcat - Vertical concatenation
function result = vertcat(this, varargin)
    % If one of the elements to be concatenated is sparse, the result will be
    % given in a sparse format (for this we'll call the sgem version of this
    % function)
    needSparse = 0;
    for i = 1:length(varargin)
        if issparse(varargin{i})
            needSparse = 1;
        end
    end
    if needSparse == 1
        result = vertcat(sparse(this), varargin{1});
        return;
    end

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
        error('Incompatible sizes for gem:vertcat');
    end

    % We make sure both objects are gem objects
    if ~isequal(class(this), 'gem') && ~isequal(class(this), 'sgem')
        this = gemify(this);
    end
    if ~isequal(class(varargin{1}), 'gem') && ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = gemify(varargin{1});
    end

    %% Now we can concatenate the two objects
    if isequal(class(this), 'gem') && isequal(class(varargin{1}), 'gem')
        newObjectIdentifier = gem_mex('vertcat', this.objectIdentifier, varargin{1}.objectIdentifier);
        result = gem('encapsulate', newObjectIdentifier);
    else
        newObjectIdentifier = sgem_mex('vertcat', this.objectIdentifier, varargin{1}.objectIdentifier);
        result = sgem('encapsulate', newObjectIdentifier);
    end
end

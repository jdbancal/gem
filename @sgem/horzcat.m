% horzcat - Horizontal concatenation
function result = horzcat(this, varargin)
    % If the concatenation involves more than 2 objects, we deal with each pair
    % independently
    if length(varargin) > 1
        result = horzcat(this, varargin{1});
        for i = 2:length(varargin)
            result = horzcat(result, varargin{i});
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
    if s1(1) ~= s2(1)
        error('Incompatible sizes for sgem:horzcat');
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
    % It seems to be quicker to perform the concatenation by this way

    % Then both matrices should be sparse already
    if ~isequal(class(this), 'sgem')
        this = sparse(this);
    end
    if ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = sparse(varargin{1});
    end

    % We resize both matrices
    this = this*sgem(1:s1(2), 1:s1(2), 1, s1(2), s1(2)+s2(2));
    varargin{1} = varargin{1}*sgem(1:s2(2), s1(2)+[1:s2(2)], 1, s2(2), s1(2)+s2(2));

    result = this + varargin{1};
    return;
    
%     % Otherwise we just use a direct method
%     if isequal(class(this), 'gem') && isequal(class(varargin{1}), 'gem')
%         newObjectIdentifier = gem_mex('horzcat', this.objectIdentifier, varargin{1}.objectIdentifier);
%         result = gem('encapsulate', newObjectIdentifier);
%     else
%         newObjectIdentifier = sgem_mex('horzcat', this.objectIdentifier, varargin{1}.objectIdentifier);
%         result = sgem('encapsulate', newObjectIdentifier);
%     end
end

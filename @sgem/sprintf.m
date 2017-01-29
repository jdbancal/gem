% sprintf - redirects to matlab's sprintf function with double conversion
function result = sprintf(varargin)

    % We convert gem objects to double
    for i = 1:nargin
        if isequal(class(varargin{i}), 'gem') || isequal(class(varargin{i}), 'sgem')
            varargin{i} = double(varargin{i});
        end
    end

    % Now we call matlab's plot function
    result = sprintf(varargin{:});
    
    if nargout == 0
        clear result;
    end
end

% sprintf - redirects to matlab's sprintf function with double conversion
function result = sprintf(varargin)

    % We convert sgem objects to dense objects, so that the gem method is
    % called instead
    for i = 1:nargin
        if isequal(class(varargin{i}), 'sgem')
            varargin{i} = full(varargin{i});
        end
    end

    % Now we call matlab's plot function
    result = sprintf(varargin{:});
    
    if nargout == 0
        clear result;
    end
end

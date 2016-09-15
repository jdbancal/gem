% plot - redirects to matlab's plot function with double conversion
function result = plot(varargin)

    % We convert gem objects to double
    for i = 1:nargin
        if isequal(class(varargin{i}), 'gem') || isequal(class(varargin{i}), 'sgem')
            varargin{i} = double(varargin{i});
        end
    end

    % Now we call matlab's plot function
    result = plot(varargin{:});
    
    if nargout == 0
        clear result;
    end
end

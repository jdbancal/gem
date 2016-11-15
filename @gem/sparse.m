% sparse - converts the matrix into sparse format
function result = sparse(varargin)
    if nargin == 1
        % We call the sgem constructor with the current object
        result = sgem(varargin{1});
        return;
    end
    
    % If we are called with several parameters, this means that we wish to
    % construct a new sgem object. So we transfer these parameters to the
    % sgem constructor.
    if nargin == 2
        result = sgem([], [], [], varargin{1}, varargin{2});
    elseif nargin == 3
        result = sgem(varargin{1}, varargin{2}, varargin{3});
    elseif nargin == 5
        result = sgem(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    else
        error('Wrong number of arguments in sgem::sparse');
    end
end

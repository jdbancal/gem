% cat - Concatenates matrices
function result = cat(dim, varargin)

switch dim
    case 1
        result = vertcat(varargin{1}, varargin{2:end});
    case 2
        result = horzcat(varargin{1}, varargin{2:end});
    otherwise
        error('Wrong dimension in sgem::cat');
end

end

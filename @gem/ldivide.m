% ldivide - element-wise left division A.\B
function result = ldivide(this, varargin)
    if length(varargin) ~= 1
        error('Wrong number of arguments in gem::ldivide');
    end

    % If the second type is more elaborated than just a number, we let the
    % corresponding class take care of performing the operation
    if ~isnumeric(varargin{1})
        result = rdivide(varargin{1}, this);
        return;
    end
    
    % We call rdivide
    result = varargin{1}./this;
end

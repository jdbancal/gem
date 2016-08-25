% ldivide - element-wise left division A.\B
function result = ldivide(this, varargin)
    if length(varargin) ~= 1
        error('Wrong number of arguments in gem::ldivide');
    end

    % We call rdivide
    result = varargin{1}./this;
end

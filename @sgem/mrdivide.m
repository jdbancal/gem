% mrdivide - right division A/B (only supported if B is a scalar
%            at the moment)
function result = mrdivide(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in sgem::mrdivide');
    end
    
    % We check that the denominator is a scalar
    if numel(varargin{1}) > 1
        error('Only right division by a scalar is supported at the moment');
    end

    result = this./varargin{1};
end

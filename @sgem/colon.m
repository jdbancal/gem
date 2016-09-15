% colon - 
%   a:b  is [a, a+1, ..., m]          with m as close to b as possible
%   a:step:b  is [a, a+step, ..., m]  with m as close to b as possible
function result = colon(varargin)
    if nargin < 2 && nargin > 3
        error('Wrong number of arguments in sgem::colon');
    end

    % We make sure no elements is sparse
    for i = 1:nargin
        % We only care about the first element of each argument
        if numel(varargin{i}) > 1
            % For some reason, matlab doesn't interpret a=a(1) as a call to
            % subsref... so we do it by hand.
            sub.type='()';
            sub.subs={[1]};
            varargin{i} = subsref(varargin{i}, sub);
        end
        varargin{i} = gem(varargin{i});
    end
    
    % We call the colon method for full gem objects
    result = colon(varargin{:});
end

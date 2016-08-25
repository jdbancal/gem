% gemRand - Uniformly distributed high precision pseudorandom numbers
%
% Usage :
%   gemRand         returns one pseudo-random number
%   gemRand(N)      returns an NxN matrix of pseudo-random numbers
%   gemRand(N,M)    returns an NxM matrix of pseudo-random numbers
%   gemRand([N,M])  idem
%
% Note : The values produced by this function depend on the seed (which can
% be set by the function gemRng), as well as on the gem working precision.
function result = gemRand(varargin)
    % Let's check the input's validity
    if nargin == 0
        size = [1, 1];
    elseif nargin == 1
        if numel(varargin{1}) == 1
            size = varargin{1}*[1 1];
        elseif numel(varargin{1}) == 2
            size = [varargin{1}(1) varargin{1}(2)];
        else
            error('Unexpected argument in gemRand');
        end
    elseif nargin == 2
        if (numel(varargin{1}) ~= 1) || (numel(varargin{2}) ~= 1) || ~isnumeric(varargin{1}) || ~isnumeric(varargin{2})
            error('Unexpected arguments in gemRand');
        else
            size = [round(varargin{1}), round(varargin{2})];
        end
    else
        error('Unexpected arguments in gemRand');
    end

    % We make sure the default working precision has been set
    gemWorkingPrecision;
    
    % We call the rand procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('rand', size);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

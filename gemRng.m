% gemRng - Setting the seed of the pseudorandom number generator gemRand
%
% Usage :
%   gemRng        uses a random seed
%   gemRng(seed)  uses seed 1+round(abs(seed))
%
% Note : The values produced by the gemRand function depend on the seed 
% used here, as well as on the gem working precision.
function gemRng(seed)
    % We check the input
    if nargin == 0
        seed = randi(2^32,1,1,'uint32');
    elseif ~isequal(class(seed),'double')
        seed = double(seed);
    end

    if (nargin > 1) || ~isnumeric(seed) || (numel(seed) ~= 1)
        error('Unexpected argument in gemRng')
    end
    
    % We call the initialization procedure
    gem_mex('randInit', 1+round(abs(seed)));
end

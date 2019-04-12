% gemRandn - high precision pseudorandom numbers distributed according to
% the standard normal distribution
%
% Usage :
%   gemRandn         returns one pseudo-random number
%   gemRandn(N)      returns an NxN matrix of pseudo-random numbers
%   gemRandn(N,M)    returns an NxM matrix of pseudo-random numbers
%   gemRandn([N,M])  idem
%
% Note : The values produced by this function depend on the seed (which can
% be set by the function gemRng), as well as on the gem working precision.
function result = gemRandn(varargin)
    % Let's check the input's validity
    if nargin == 0
        size = [1, 1];
    elseif nargin == 1
        if numel(varargin{1}) == 1
            size = varargin{1}*[1 1];
        elseif numel(varargin{1}) == 2
            size = [varargin{1}(1) varargin{1}(2)];
        else
            error('Unexpected argument in gemRandn');
        end
    elseif nargin == 2
        if (numel(varargin{1}) ~= 1) || (numel(varargin{2}) ~= 1) || ~isnumeric(varargin{1}) || ~isnumeric(varargin{2})
            error('Unexpected arguments in gemRandn');
        else
            size = [round(varargin{1}), round(varargin{2})];
        end
    else
        error('Unexpected arguments in gemRandn');
    end
    nbElements = prod(size);
    
    % We use Marsaglia's polar method (adapted and corrected from https://www.projectrhea.org/rhea/index.php/The_principle_for_how_to_generate_random_samples_from_a_Gaussian_distribution )

    % We first get uniformly random numbers
    sel = [];
    while (length(sel)*2 < nbElements)
        if (~isempty(sel))
            disp('trying again');
        end
        w1 = 2*gemRand(1, ceil(nbElements*0.7)) - 1;
        w2 = 2*gemRand(1, ceil(nbElements*0.7)) - 1;
        r = w1.^2 + w2.^2;
        sel = find(r < 1);
    end
    
    % We only use the number minimum number of elements
    sel = sel(1:ceil(nbElements/2));
    w1 = w1(sel);
    w2 = w2(sel);
    r = r(sel);
    
    % Now we obtain the random numbers
    r = sqrt(-2*log(r)./r);
    result = [w1.*r, w2.*r];
    
    % And format them in the desired way
    result = reshape(result(1:nbElements), size);
end

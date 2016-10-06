% precision = gemDisplayPrecision(precision)
% 
% Getting and setting the display precision for all gem objects.
%
% To get the current precision used when displaying gem objects, use
%   precision = gemDisplayPrecision
% To display all gem objects with 30 digits, use
%   gemDisplayPrecision(30)
function precision = gemDisplayPrecision(precision)
    tmp = gem(0);
    if nargin < 1
        precision = tmp.getDisplayPrecision;
    else
        tmp.setDisplayPrecision(precision);
        % We also set the precision for the sparse objects
        stmp = sgem(0);
        stmp.setDisplayPrecision(precision);
    end
end

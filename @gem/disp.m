% display function, without printing of name
%
% possible usages:
%   disp(a)                  prints a with default precision
%   disp(a, precision)       prints a with given precision
%
% Note: Here, the precision is the number of printed digits (including leading zeros)
function disp(this, arg2)
    % The default display precision
    displayPrecision = this.getDisplayPrecision;
    
    % Some versions of matlab provide the name as a second argument
    if (nargin >= 2)
        displayPrecision = double(arg2);
    end
    
    gem_mex('display', this.objectIdentifier, displayPrecision);
    disp(' ');
end

% precision = gemWorkingPrecision(precision)
% 
% Getting and setting the working precision for all gem objects.
%
% To get the current working precision , use
%   precision = gemWorkingPrecision
% To set the default precision for all newly created gem object to 100, use
%   gemWorkingPrecision(100)
%
% The working precision is the precision at which new gem object are 
% created. This precision applies to objects created from matlab, but also 
% to all result of operations performed by the gem library. Therefore, you 
% should always start your work by calling this function with the level of
% precision with which you want computations to be performed. Failure to
% call this function will result in a default precision of 50 digits being
% used by the library.
%
% Example : calling gemWorkingPrecision(10) yields the following behavior:
%   gem(3.1,3)*gem(2,3)  ->  6.203125000
% Here, additional digits have been added to the result.
% These digits are zeros if both numbers in the product are integers:
%   gem(312,3)*gem(223,3)  ->  69576
% Zeros after the decimal points are not printed for integer numbers, so
% we can't see them here. To make sure there are 10 digits in total, we
% can display it with a fractional number :
%   [1.1 gem(312,3)*gem(223,3)] ->  1.100000000  69576.00000
% this shows that the result of the multiplication finishes with exactly
% 5 zeros.
%
% The working precision also puts an upper bound on the number of digits 
% kept after an operation. Therefore
%   gem(1.23456789012345,15)*gem(1.23456789012345,15)
% only returns 10 digits : 1.524157875
function precision = gemWorkingPrecision(precision)
    tmp = gem(0);
    if nargin < 1
        precision = tmp.getWorkingPrecision;
    else
        tmp.setWorkingPrecision(precision);
        % We also set the precision for the sparse objects
        stmp = sgem(0);
        stmp.setWorkingPrecision(precision);
    end
end

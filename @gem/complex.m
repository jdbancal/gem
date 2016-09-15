% complex(a,b) - Returns a + 1i*b
% 
% Creates a complex matrix from real and imaginary components a and b
function result = complex(a, b)
    % If no second argument is given, the output is the real part
    if length(varargin) ~= 1
	result = a;
        return;
    end

    result = a + 1i*b;
end

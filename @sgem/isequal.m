% isequal - True if arrays are numerically equal (both in type, size,
%           digits and in precision)
%           Sparsity is not required to match.
%
% example : isequal(sgem([1 0 1.2]), [1 0 1.2]) gives 1
function result = isequal(varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) < 2
        error('Wrong number of arguments in sgem::isequal');
    end
    
    % If there are more than two items to compare, we iterate the
    % comparison on successive pairs
    if length(varargin) > 2
        result = zeros(1,length(varargin)-1);
        for i = 1:length(varargin)-1
            result(i) = isequal(varargin{i}, varargin{i+1});
        end
        result = (sum(result) == length(varargin)-1);
        return;
    end
    % From now one we only deal with two objects
    
    % First, we check if the two objects have the same dimension
    if ~isequal(size(varargin{1}), size(varargin{2}))
        result = false;
        return;
    end
    
    % Next, we check whether both objects are either real or complex
    if isreal(varargin{1}) + isreal(varargin{2}) == 1
        result = false;
        return;
    end

% We allow sparse object to be equal to full ones (so no check of sparsity)
%     if ~issparse(varargin{1}) || ~issparse(varargin{2})
%         result = false;
%         return;
%     end

    % Now, we make sure that both objects are sgem objects
    if ~isequal(class(varargin{1}), 'sgem')
        varargin{1} = sgem(varargin{1});
    elseif ~isequal(class(varargin{2}), 'sgem')
        varargin{2} = sgem(varargin{2});
    end
    
    % Now we check whether the precision of both matrices match
    if ~isequal(precision(varargin{1}), precision(varargin{2}))
        result = false;
        return;
    end
    
    % Finally, we check whether the numerical values match
    result = logical(sgem_mex('identicalValues', varargin{1}.objectIdentifier, varargin{2}.objectIdentifier));
end

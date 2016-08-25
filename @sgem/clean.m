% clean - removes small real and imaginary components
%
% Supported formats:
%   y = clean(x, tol)  Removes all components which are smaller than tol
%                      in absolute value
%   y = clean(x)       Same as above with tol = 10^(-gemWorkingPrecision)
function result = clean(this, tol)
    if nargin < 2
        % The default tolerance
        tol = 10^(-gemWorkingPrecision);
    end
    
    % We make sure tol is a single gem number
    if isequal(class(tol), 'sgem')
        tol = gem(tol);
    end
    if ~isequal(class(tol), 'gem')
        if ~isnumeric(tol)
            error('The tolerance should be a number in sgem::clean');
        end
        tol = gem(tol);
    end
    if numel(tol) ~= 1
        error('The tolerance should be a single number in sgem::clean');
    end
    if tol < 0
        error('The tolerance should be a positive number in sgem::clean');
    end

    % We call the cleaning procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = sgem_mex('clean', this.objectIdentifier, tol.objectIdentifier);

    % ...  and create a new matlab object to keep this handle
    result = sgem('encapsulate', newObjectIdentifier);
end

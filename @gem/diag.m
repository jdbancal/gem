% diag - returns the diagonal elements of a matrix
%
% supported formats:
%   diag(M) : the main diagonal of matrix M
%   diag(M,k) : the k-th diagonal of matrix M (k=0 is the main diagonal
%   diag(d) : a matrix with diagonal d
%   diag(d,k) : a matrix with k-th diagonal d
function result = diag(this, k)
    if (nargin < 2)
        k = 0;
    else
        k = double(k);
        if (numel(k) ~= 1) || (round(k) ~= k)
            error('Unexpected argument in gem::diag');
        end
    end

    % diag of empty is empty
    if isempty(this)
        result = this;
        return;
    end
    
    % We check that the user is not trying to create a matrix
    s = size(this);
    if min(s) == 1
        % We want the indices in a column
        if s(2) > s(1)
            this = this.';
        end
        % We call the angle procedure. Since the function creates a
        % new object with the result, we keep the corresponding handle...
        newObjectIdentifier = gem_mex('diagCreate', this.objectIdentifier, k);
    else
        if (k+1 > s(2)) || (-k+1 > s(1))
            error('diagonal index is larger than the matrix');
        end

        % We call the angle procedure. Since the function creates a
        % new object with the result, we keep the corresponding handle...
        newObjectIdentifier = gem_mex('diagExtract', this.objectIdentifier, k);
    end
    

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end

% find - Find indices of nonzero elements
%
% Supported formats:
%   I = find(x)            Returns the linear indices I for which x(I) ~= 0
%   I = find(x,k)          Returns the linear indices of at most k of the
%                          first non-zero elements in x
%   I = find(x,k,'first')  Same as above
%   I = find(x,k,'last')   Rather returns at most k of the last nonzero
%                          elements in x
%   I = find(x)            Returns the linear indices I for which x(I) ~= 0
%   [I, J] = find(x)       Returns the indices for which x(I(k),J(k)) ~= 0
%   [I, J, V] = find(x)    Also returns the nonzero values V
function [I, J, V] = find(this, k, mode)
    if (nargin >= 2) && ((numel(k) ~= 1) || (~isnumeric(k)))
        error('Unknown parameter in sgem::find');
    end
    if (nargin >= 3) && (~isequal(lower(mode), 'first')) && (~isequal(ower(mode), 'last'))
        error('Unknown parameter in sgem::find');
    end

    [I, J, V] = sgem_mex('find', this.objectIdentifier);
    % Indices in matlab start from 1
    I = I+1;
    J = J+1;
    % The third result is a gem object
    V = gem('encapsulate',V);

    if (nargin >= 2)
        k = min(k, length(I));
        if (nargin < 3) || (isequal(lower(mode), 'first'))
            I = I(1:k);
            J = J(1:k);
            V = V(1:k);
        else
            I = I(end-k:end);
            J = J(end-k:end);
            V = V(end-k:end);
        end
    end

    % For row vectors, the outputs are row vectors
    if size(this,1) == 1
        I = I';
        J = J';
        V = V';
    end
    if nargout <= 1
        % We translate the indices into linear indices
        I = size(this,1)*(J-1) + I;
    end
end

% eig - eigenvalues and eigenvectors
%
% supported formats :
%   e = eig(a) : eigenvalues of a
function [V D] = eig(this, varargin)
    % This function can involve up to two argument
    if length(varargin) > 1
        error('Wrong number of arguments in gem::eig');
    end
    
    % We check how many outputs are 
    if nargout <= 2
        [newObjectIdentifierD newObjectIdentifierV] = gem_mex('eig', this.objectIdentifier);

        % ...  and create a new matlab object to keep this handle
        D = gem('encapsulate', newObjectIdentifierD);
        V = gem('encapsulate', newObjectIdentifierV);
        
        % Now we check if D has a block structure
        blocks = (1+find(diag(D,1)))/2;
        
        % Eigenvalues in blocks are doubled

        % The other ones are simple
        
        if nargout <= 1
            V = diag(D);
        end
    else
        error('Unsupported call to gem::eig')
    end

end

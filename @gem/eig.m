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
        
        % We normalize the eigenvectors (this should not be done if the
        % option 'nobalance' is passed (once this option is implemented).
        V = V*diag(1./sqrt(diag(V'*V)));

        % To check whether the eigendecomposition is correct:
        % V*D*inv(D)-(this)
%         precision = double(abs(norm(mtimes(V,mtimes(D,inv(V))) - this,1)));
%         disp(['Eigenvalue decomposition precision: ', num2str(precision)]);
%         if (precision > 1e-30)
%             disp('WARNING : big eigendecomposition error!!!');
%         end
        
        if nargout <= 1
            V = diag(D);
        end
        
    else
        error('Unsupported call to gem::eig')
    end

end

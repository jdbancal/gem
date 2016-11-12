% eigs - partial eigenvalues and eigenvectors
%
% supported formats :
%   e = eigs(a, [k])          : k eigenvalues of with largest magnitude
%   [v d] = eigs(a, [k])      : k eigenvectors and eigenvalues of a 
%                               with largest magnitude (a*v = v*d)
%   [v d] = eigs(a, k, sigma) : first k eigenvectors and eigenvalues of a
%                               with sigma = 'lm', 'sm' for the eigenvalues
%                               with largest or smallest magnitude
% 
% The default value of k is 6.
%
% Note:
%  - The precision of the computation performed here is gemWorkingPrecision*2/3
%  - It appears that sometimes the order of two eigenvalues can be
%    inverted, so if this is crucial information, one should double check 
%    the last eigenvalues by computing a few more eigenvalues than needed.
function [V D] = eigs(this, varargin)
    % This function can involve at most two parameters
    if length(varargin) > 2
        error('Wrong number of arguments in sgem::eigs');
    end
    
    if (length(varargin) > 0) && (~isnumeric(varargin{1}) || (numel(varargin{1}) ~= 1))
        error('The second argument of sgem::eigs must be a single number');
    elseif (length(varargin) > 0) && ~isequal(class(varargin{1}), 'double')
        % We make sure that the number of eigenvalues to be computed was
        % specified as a double
        varargin{1} = double(varargin{1});
    end
    
    % Extract the requested number of eigenvalues
    if length(varargin) > 0
        nbEigenvalues = varargin{1};
    else
        nbEigenvalues = min(size(this,1), 6);
    end
    
    % The number of eigenvalues computed must be larger than zero
    if nbEigenvalues < 1
        error('sgem::eigs cannot compute less than 1 eigenvalue');
    end
    
    if nbEigenvalues > size(this,1) - 2 + ishermitian(this)
        error('Too many eigenvalues for sgem::eigs');
    end
    
    % We check if there is a second parameter
    if (length(varargin) > 1) && (~ischar(varargin{2}))
        error('The third argument of sgem::eigs must be a text');
    end

    % We extract the requested type of eigenvalues
    type = 1;
    if length(varargin) > 1
        switch lower(varargin{2})
            case 'lm'
                % Largest magnitude
                type = 1;
            case 'sm'
                % Smallest magnitude
                type = 2;
% % The following cases are not implemented at the moment:
%             case 'lr'
%                 % Largest real component
%             case 'sr'
%                 % Smallest real component
%             case 'li'
%                 % Largest imaginary component
%             case 'si'
%                 % Smallest imaginary component
            otherwise
                error('Third argument of sgem::eigs not recognized');
        end
    end
    
    % We check how many outputs are
    if nargout <= 2
        % The matrix must be square
        if size(this, 1) ~= size(this,2)
            error('Matrix must be square in sgem::eigs');
        end
        
        [newObjectIdentifierV newObjectIdentifierD] = sgem_mex('eigs', this.objectIdentifier, nbEigenvalues, type);

        % ...  and create a new matlab object to keep this handle
        V = gem('encapsulate', newObjectIdentifierV);
        D = gem('encapsulate', newObjectIdentifierD);
        
        % We normalize the eigenvectors (this should not be done if the
        % option 'nobalance' is passed (once this option is implemented).
        V = V*diag(1./sqrt(diag(V'*V)));

        if nargout <= 1
            V = diag(D);
        end
        
    else
        error('Unsupported call to sgem::eigs')
    end

end

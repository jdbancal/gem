% eigs - partial eigenvalues and eigenvectors
%
% supported formats :
%   e = eigs(a, [k])          : k eigenvalues of with largest magnitude
%   [v d] = eigs(a, [k])      : k eigenvectors and eigenvalues of a 
%                               with largest magnitude (a*v = v*d)
%   [v d] = eigs(a, k, sigma) : first k eigenvectors and eigenvalues of a
%                               smallest in magnitude to sigma.
%                               With sigma = 'lm', 'sm', eigenvalues
%                               with largest or smallest magnitude.
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
        error('Wrong number of arguments in gem::eigs');
    end
    
    if (length(varargin) > 0) && (~isnumeric(varargin{1}) || (numel(varargin{1}) ~= 1))
        error('The second argument of gem::eigs must be a single number');
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
        error('gem::eigs cannot compute less than 1 eigenvalue');
    end
    
    % We check if there is a second parameter
    type = 1;
    sigma = 0;
    if length(varargin) > 1
        if ischar(varargin{2})
            switch lower(varargin{2})
                case 'lm'
                    % Largest magnitude
                    type = 1;
                case 'sm'
                    % Smallest magnitude
                    type = 2;
% % The following cases are not implemented at the moment:
%                 case 'lr'
%                     % Largest real component
%                 case 'sr'
%                     % Smallest real component
%                 case 'li'
%                     % Largest imaginary component
%                 case 'si'
%                     % Smallest imaginary component
                otherwise
                    error('Third argument of gem::eigs not recognized');
            end
        else
            if numel(varargin{2}) > 1
                error('sigma must be a single number in gem::eigs');
            end
            type = 2;
            sigma = varargin{2};
        end
    end
    
    % There should be no more parameters
    if length(varargin) > 2
        error('Too many arguments in gem::eigs');
    end

    if nbEigenvalues > size(this,1) - 2 + ishermitian(this)
        if nbEigenvalues > size(this,1)
            error('Too many eigenvalues asked for in gem::eigs');
        end
                
        % We use eig to compute all eigenvalues
%        warning('Too many eigenvalues for eigs, using eig instead.');
        if nargout == 2
            [V D] = eig(this);
            % We make sure the eigenvalues are ordered (not guaranteed by
            % eig)
            [D reorder] = sort(diag(D),1,'descend');
            D = diag(D);
            subV.type='()';
            subV.subs={':' reorder};
            V = subsref(V, subV);
            if isequal(type, 2)
                subV.type='()';
                subV.subs={[1:size(V,1)] [size(V,2):-1:size(V,2)-nbEigenvalues+1]};
                V = subsref(V, subV);
                subD.type='()';
                subD.subs={[size(D,1):-1:size(D,1)-nbEigenvalues+1] [size(D,2):-1:size(D,2)-nbEigenvalues+1]};
                D = subsref(D, subD);
            elseif nbEigenvalues < size(D,1)
                subV.type='()';
                subV.subs={[1:size(V,1)] [1:nbEigenvalues]};
                V = subsref(V, subV);
                subD.type='()';
                subD.subs={[1:nbEigenvalues] [1:nbEigenvalues]};
                D = subsref(D, subD);
            end
        else
            V = eig(this);
            % We make sure the eigenvalues are ordered (not guaranteed by
            % eig)
            V = sort(V,1,'descend');
            if isequal(type, 2)
                subV.type='()';
                subV.subs={[size(V,1):-1:size(V,1)-nbEigenvalues+1] [1]};
                V = subsref(V, subV);
            elseif nbEigenvalues < size(V,1)
                subV.type='()';
                subV.subs={[1:nbEigenvalues] [1]};
                V = subsref(V, subV);
            end
        end
        return;
    end
    
    % We check how many outputs are
    if nargout <= 2
        % The matrix must be square
        if size(this, 1) ~= size(this,2)
            error('Matrix must be square in gem::eigs');
        end
        
        % We make sure sigma is a gem object
        if ~isequal(class(sigma),'gem')
            sigma = gem(sigma);
        end

        % The algorithm we use doen't handle eigenvalues equal to sigma. So
        % this variable tells how many times it needs to be added
        % mannually.
        additionalSigmaMultiplicity = 0;
        
        % We make sure we won't try to find non-zero eigenvalues when there
        % are no more (this is numerically unstable)
        if isequal(type,1)
            rankMatrix = rank(this);
            if rankMatrix < nbEigenvalues
                if nargout >= 2
                    error('Eigenvectors are not computed for null eigenvalues.')
                end
                if rankMatrix == 1
                    warning('There is only one non-zero eigenvalues, computing this one only.');
                else
                    warning(['There are only ', num2str(rankMatrix), ' non-zero eigenvalues, computing these ones only.']);
                end
                additionalSigmaMultiplicity = nbEigenvalues-rankMatrix;
                nbEigenvalues = rankMatrix;
            end
        end
        
        % We make sure we won't try to invert a singular matrix (this is
        % numerically unstable as well)
        if isequal(type,2)
            rankShifted = rank(this-sigma*eye(size(this)));
            if size(this,1) - rankShifted >= nbEigenvalues
                % Then we know that all requested eigenvalues are equal to
                % sigma, but we still haven't computed corresponding
                % eigenvectors (these would be given by a function such as
                % 'null')
                V = sigma*ones(nbEigenvalues,1);
                if nargout >= 2
                    error('Eigenvectors are not computed when sigma is an eigenvalue.')
                end
                return;
            elseif size(this,1) - rankShifted > 0
                error('Sigma is an eigenvalue of the considered matrix. Consider perturbing it a little bit to allow the numerical method to run smoothly.')
            end
        end
        
        [newObjectIdentifierV newObjectIdentifierD] = gem_mex('eigs', this.objectIdentifier, nbEigenvalues, type, sigma.objectIdentifier);

        % ...  and create a new matlab object to keep this handle
        V = gem('encapsulate', newObjectIdentifierV);
        D = gem('encapsulate', newObjectIdentifierD);
        
        % We normalize the eigenvectors (this should not be done if the
        % option 'nobalance' is passed (once this option is implemented)).
        V = V*diag(1./sqrt(diag(V'*V)));
        
        if nargout <= 1
            V = diag(D);
            if additionalSigmaMultiplicity > 0
                V = [V; sigma*ones(additionalSigmaMultiplicity,1)];
            end
        end
        
    else
        error('Unsupported call to gem::eigs')
    end

end

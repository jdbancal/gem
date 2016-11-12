% svds - Partial singular value decomposition
%
% supported formats :
%   s = svds(a, [k])            : k largest singular values of a
%   [u s v] = svds(a, [k])      : k largest singular values of a with left
%                                 and right economic singular vectors u and
%                                 v: s = u'*a*v. The columns of u and v are 
%                                 orthonormal, i.e. u'*u = id, v'*v = id.
%   [u s v] = svds(a, k, sigma) : sigma = 'largest', 'smallest' for the
%                                 largest or smallest singular values
%
% The default value of k is 6.
function [U S V] = svds(this, varargin)
    % This function can involve at most two parameters
    if length(varargin) > 2
        error('Wrong number of arguments in sgem::svds');
    end
    
    if (length(varargin) > 0) && (~isnumeric(varargin{1}) || (numel(varargin{1}) ~= 1))
        error('The second argument of sgem::svds must be a single number');
    elseif (length(varargin) > 0) && ~isequal(class(varargin{1}), 'double')
        % We make sure that the number of eigenvalues to be computed was
        % specified as a double
        varargin{1} = double(varargin{1});
    end
    
    % Extract the requested number of singular values
    if length(varargin) > 0
        nbSingularvalues = varargin{1};
    else
        nbSingularvalues = min(size(this,1), 6);
    end
    
    % The number of singular values computed must be larger than zero
    if nbSingularvalues < 1
        error('sgem::svds cannot compute less than 1 singular value');
    end
    
    if nbSingularvalues > size(this,1) - 2 + ishermitian(this)
        % We use svd to compute all singular values
        error('Too many singular values for sgem::svds.');
    end
    
    % We check if there is a second parameter
    if (length(varargin) > 1) && (~ischar(varargin{2}))
        error('The third argument of sgem::svds must be a text');
    end

    % We extract the requested type of singular values
    type = 'lm';
    if length(varargin) > 1
        switch lower(varargin{2})
            case 'largest'
                % Largest singular values
                type = 'lm';
            case 'smallest'
                % Smallest singular values
                type = 'sm';
            otherwise
                error('Third argument of sgem::svds not recognized');
        end
    end
    
    % We perform the computation
    if nargout <= 1
        % We only compute the eigenvalues of a*a'
        vals = eigs(this*this', nbSingularvalues, type);
        U = sqrt(vals);

        % We make sure the order of the singular values is decreasing
        if isequal(type,'sm')
            subU.type='()';
            subU.subs={[size(U,1):-1:1]};
            U = subsref(U, subU);
        end
    elseif nargout <= 2
        % We compute the eigenvalues and eigenvectors of a*a'
        [U valsU] = eigs(this*this', nbSingularvalues, type);
        S = sqrt(valsU);

        % We make sure the order of the singular values is decreasing
        if isequal(type,'sm')
            subU.type='()';
            subU.subs={[1:size(U,1)] [size(U,2):-1:1]};
            U = subsref(U, subU);
            subS.type='()';
            subS.subs={[size(S,1):-1:1] [size(S,2):-1:1]};
            S = subsref(S, subS);
        end
    elseif nargout <= 3
        % We compute the eigenvalues on both sides
        [U valsU] = eigs(this*this', nbSingularvalues, type);
        [V valsV] = eigs(this'*this, nbSingularvalues, type);
        S = sqrt(valsU);
        
        % We check that the same singular values were found on both sides
        if max(abs(valsU-valsV)) > 10^(-gemWorkingPrecision*2/3 + 2)
            warning('The left and right singular values don''t match in sgem::svds')
        end
        
        % We correct the phases of V
        V = (V.*(ones(size(V,1),1)*exp(-1i*diag(angle(U'*this*V))).'));

        % We make sure the order of the singular values is decreasing
        if isequal(type,'sm')
            subU.type='()';
            subU.subs={[1:size(U,1)] [size(U,2):-1:1]};
            U = subsref(U, subU);
            subS.type='()';
            subS.subs={[size(S,1):-1:1] [size(S,2):-1:1]};
            S = subsref(S, subS);
            subV.type='()';
            subV.subs={[1:size(V,1)] [size(V,2):-1:1]};
            V = subsref(V, subV);
        end
    else
        error('Unsupported call to sgem::svds')
    end

end

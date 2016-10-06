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

        if nargout == 1
            V = diag(D);
        end
        
%         % The following has been now coded into the c++ library
%         % We extract the complex information
%         if isreal(this)
%             % Now we check if D has a block structure
%             blocks = find(diag(D,1));
% 
%             % Eigenvalues in 1x1 blocks are simple
%             simple = setdiff(1:3,[blocks blocks+1]);
% 
%             % Eigenvalues in 2x2 blocks are doubled
%             sub1.type = '()';
%             sub1.subs = {':',blocks};
%             sub2.type = '()';
%             sub2.subs = {':',blocks+1};
%             V = subsasgn(V, sub1, subsref(V,sub1) + 1i*subsref(V,sub2));
%             V = subsasgn(V, sub2, conj(subsref(V,sub1)));
% 
%             % The other ones are simple
%             if nargout == 1
%                 V = diag(D);
%             end
%             
% 
%         else
%             % For complex matrices, we extract one eigenvector and one
%             % eigenvalue per block.
% 
%             % We normalize the eigenvectors (this should not be done if the
%             % option 'nobalance' is passed (once this option is implemented).
%             %V = V(1:end/2,1:2:end)-1i*V(end/2+1:end,1:2:end);
%             V = V*diag(1./sqrt(diag(V'*V)));
%             
%             % If we want to check that the decomposition is correct
%             %   vv*dd*inv(vv)-m
%             
%             if nargout == 1
%                 V = diag(D);
%             end
%         end
        
    else
        error('Unsupported call to gem::eig')
    end

end

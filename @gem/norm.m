% norm - computes the matrix norm
% 
% norm(x) is the 2-norm of x
function result = norm(this, p)
    % We define the kind of norm we are interested in if this was not
    % specified
    if nargin < 2
        p = 2;
    end
    
    % Some norms are only supported for vectors
    if min(size(this)) == 1
        if isequal(p,'fro') || ~isnumeric(p) || (p <= 0 && p ~= -Inf)
            error('gem::norm : unsupported vector norm');
        end
    else
        if (~isequal(p, 1)) && (~isequal(p, Inf)) && (~isequal(p,'fro')) % && (~isequal(p, 2)) % to be supported soon
            error('gem::norm : unsupported matrix norm');
        end
    end
    
    % We computing the requires norm, possibly not in the most efficient
    % way
    if min(size(this)) == 1
        if (p == 2) || (p == gem(2))
            result = sqrt(sum(this.^gem(2)));
        elseif (p == Inf) || (p == gem(Inf))
            result = max(abs(this));
        elseif (p == -Inf) || (p == gem(-Inf))
            result = min(abs(this));
        else
            result = sum(abs(this).^gem(p)).^(1/gem(p));
        end
    else
        switch p
            case 1
                result = max(sum(this));
            case Inf
                result = max(sum(this,2));
            case 'fro'
                result = sum(sum(abs(this).^gem(2)))^(1/gem(2));
        end
    end
end

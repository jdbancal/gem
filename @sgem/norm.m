% norm - computes the matrix norm
% 
% supported formats :
%   norm(x)    : norm(x, 2)
%   norm(x, p) : with p = 1, 2, -Inf, Inf for vector x
%   norm(x, p) : with p = 1, 2, Inf for matrix x
function result = norm(this, p)
    % We define the kind of norm we are interested in if this was not
    % specified
    if nargin < 2
        p = 2;
    end
    
    % Some norms are only supported for vectors
    if min(size(this)) == 1
        if isequal(p,'fro') || ~isnumeric(p) || (p <= 0 && p ~= -Inf)
            error('sgem::norm : unsupported vector norm');
        end
    else
        if (~isequal(p, 1)) && (~isequal(p, Inf)) && (~isequal(p,'fro')) && (~isequal(p, 2))
            error('sgem::norm : unsupported matrix norm');
        end
    end
    
    % We computing the requires norm, possibly not in the most efficient
    % way
    if min(size(this)) == 1
        if double(p) == 2
            result = sqrt(sum(abs(this).^gem(2)));
        elseif double(p) == Inf
            result = max(abs(this));
        elseif double(p) == -Inf
            result = min(abs(this));
        else
            result = sum(abs(this).^gem(p)).^(1/gem(p));
        end
    else
        switch p
            case 1
                result = max(sum(abs(this)));
            case 2
                result = svds(this, 1, 'largest');
            case Inf
                result = max(sum(abs(this),2));
            case 'fro'
                result = sum(sum(abs(this).^gem(2)))^(1/gem(2));
        end
    end
end

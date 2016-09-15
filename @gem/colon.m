% colon - 
%   a:b  is [a, a+1, ..., m]          with m as close to b as possible
%   a:step:b  is [a, a+step, ..., m]  with m as close to b as possible
function result = colon(varargin)
    if nargin < 2 && nargin > 3
        error('Wrong number of arguments in gem::colon');
    end

    if nargin == 2
        a = varargin{1};
        b = varargin{2};

        % In case these are matrices, we only care about their first elements
        if numel(a) > 1
            % For some reason, matlab doesn't interpret a=a(1) as a call to
            % subsref... so we do it by hand.
            sub.type='()';
            sub.subs={[1]};
            a = subsref(a, sub);
        end
        if numel(b) > 1
            sub.type='()';
            sub.subs={[1]};
            b = subsref(b, sub);
        end

        % Now we also check whether both objects are of type gem,
        % otherwise we convert them
        if ~isequal(class(a), 'gem')
            a = gem(a);
        end
        if ~isequal(class(b), 'gem')
            b = gem(b);
        end

        if a > b
            % In this case the output is empty
            result = gem([]);
            return;
        end

        % We generate the output vector
        result = [a zeros(1,double(fix(b-a)))];
        % For some reason, matlab doesn't interpret result(1)=a as a call to
        % subsasgn... so we do it by hand.
        sub.type='()';
        sub.subs={[1]};
        result = subsasgn(result, sub, a);
        for i = 1:double(fix(b-a));
            sub.subs={[1+i]};
            result = subsasgn(result, sub, a+i);
        end
    else
        a = varargin{1};
        step = varargin{2};
        b = varargin{3};

        % In case these are matrices, we only care about their first elements
        if numel(a) > 1
            % For some reason, matlab doesn't interpret a=a(1) as a call to
            % subsref... so we do it by hand.
            sub.type='()';
            sub.subs={[1]};
            a = subsref(a, sub);
        end
        if numel(step) > 1
            sub.type='()';
            sub.subs={[1]};
            step = subsref(step, sub);
        end
        if numel(b) > 1
            sub.type='()';
            sub.subs={[1]};
            b = subsref(b, sub);
        end

        % Now we also check whether both objects are of type gem,
        % otherwise we convert them
        if ~isequal(class(a), 'gem')
            a = gem(a);
        end
        if ~isequal(class(step), 'gem')
            step = gem(step);
        end
        if ~isequal(class(b), 'gem')
            b = gem(b);
        end

        if (a > b) || (step == 0)
            % In this case the output is empty
            result = gem([]);
            return;
        end

        % We generate the output vector
        result = gem(zeros(1,double(1+fix((b-a)/step))));
        % For some reason, matlab doesn't interpret result(1)=a as a call to
        % subsasgn... so we do it by hand.
        sub.type='()';
        sub.subs={[1]};
        result = subsasgn(result, sub, a);
        for i = 1:double(fix((b-a)/step));
            sub.subs={[1+i]};
            result = subsasgn(result, sub, a+i*step);
        end
    end

end

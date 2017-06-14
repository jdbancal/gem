% sprintf - redirects to matlab's sprintf function with double conversion
function result = sprintf(varargin)

    if (nargin == 2) && ischar(varargin{1}) ...
            && (varargin{1}(1) == '%') && (varargin{1}(2) ~= '#') && (sum(varargin{1}=='.')==1) && (varargin{1}(end) == 'g')

        % Then we construct the string by hand
        
        % parameters
        format = varargin{1};
        width = abs(str2double(format(2:find(format=='.')-1)));
        precision = str2double(format(find(format=='.')+1:end-1));
        width = max(width,precision+1);
        
        txt = toStrings(varargin{2},precision);
        if ~iscell(txt)
            tmp = txt;
            txt = cell(1);
            txt{1} = tmp;
        end
        for i = 1:numel(txt)
            if length(txt{i}) < width
                txt{i} = [txt{i}, char(32*ones(1,width-length(txt{i})))];
            end
        end
        
        text = '';
        for i = 1:size(txt,1)
            text(i,:) = [txt{i,:}];
        end
        
        result = text;
        
    else
        % We convert gem objects to double and rely on matlab's
        % implementation
        for i = 1:nargin
            if isequal(class(varargin{i}), 'gem') || isequal(class(varargin{i}), 'sgem')
                varargin{i} = double(varargin{i});
            end
        end

        % Now we call matlab's plot function
        result = sprintf(varargin{:});
    end

    if nargout == 0
        clear result;
    end

end

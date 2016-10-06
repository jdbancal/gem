% subsref - selects part of the matrix
function result = subsref(this, varargin)


% We do some checks
if isequal(varargin{1}, '()') && ((length(varargin) > 1) || (length(varargin{1}.subs) > 2))
    error('Wrong call to gem::subsref');
end


switch varargin{1}(1).type
    case '()'
        % We extract the coordinates of the matrix
        indices = cell(1,length(varargin{1}.subs));
        s = size(this);
        for i = 1:length(varargin{1}.subs)
            if isnumeric(varargin{1}.subs{i})
                indices{i} = varargin{1}.subs{i};
            elseif isequal(varargin{1}.subs{i},':')
                if length(varargin{1}.subs) == 1
                    % We are calling with a single index a(:)
                    indices{i} = [1:prod(s)]';
                else
                    % We are calling with several indices, as in a(:,1)
                    indices{i} = 1:s(i);
                end
            else
                error('Unrecognized indexing in gem::subsref')
            end
        end
    case '.'
        % Somehow, now that we overloaded the subsref function, methods are
        % not accessible anymore... we restore it for the methods meant to
        % be public
        switch varargin{1}(1).subs
            case 'getWorkingPrecision'
                result = this.getWorkingPrecision;
                return;
            case 'setWorkingPrecision'
                this.setWorkingPrecision(varargin{1}(2).subs{1});
                return;
            case 'getDisplayPrecision'
                result = this.getDisplayPrecision;
                return;
            case 'setDisplayPrecision'
                this.setDisplayPrecision(varargin{1}(2).subs{1});
                return;
            case 'objectIdentifier'
                % The sgem class is allowed to access this private
                [ST I] = dbstack('-completenames');
                if (length(ST) < 2) || isempty(strfind(ST(2).file,'/@sgem/'))
                    error('Only sgem.m is allowed to access this property.');
                end
                result = this.objectIdentifier;
                return;
            otherwise
                error('Unsupported referencing in gem::subsref');
        end
    otherwise
        error('Wrong call to gem::subsref');
end


%% Let's make some checks
if length(indices) == 1
    if isempty(indices{1})
        % The result is an empty matrix
        result = gem([]);
        return;
    end
    if length(size(indices{1})) > 2
        error('Indices can have at most two dimensions in gem::subsref');
    end
    if (min(min(indices{1})) < 1) || (max(max(indices{1})) > prod(s))
        error('Indices out of bound in gem::subsref');
    end
elseif length(indices) == 2
    if (min(size(indices{1})) ~= 1) || (length(size(indices{1})) > 2) || (min(size(indices{2})) ~= 1) || (length(size(indices{2})) > 2)
        error('Invalid indices in gem::subsref');
    end
    if (min(indices{1}) < 1) || (max(indices{1}) > s(1)) || (min(indices{2}) < 1) || (max(indices{2}) > s(2))
        error('Indices out of bound in gem::subsref');
    end
else
   error('Gem objects are at most two dimensional in gem::subsref');
end


%% Now we extract the requested numbers

if length(indices) == 1
    % Then indices have been specified for one dimension, as in a(1:2)
    % so we call the subsref procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...

    newObjectIdentifier = gem_mex('subsref', this.objectIdentifier, indices{1}-1);
    
    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
else
    % Then indices have been specified for both dimensions, as in a(1:2,:)
    % so we call the subsref procedure. Since the function creates a
    % new object with the result, we keep the corresponding handle...
    newObjectIdentifier = gem_mex('subsref', this.objectIdentifier, indices{1}-1, indices{2}-1);

    % ...  and create a new matlab object to keep this handle
    result = gem('encapsulate', newObjectIdentifier);
end




%result = prod(size(this));
end

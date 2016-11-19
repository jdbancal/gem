% mpower - matrix power A^B -- only supported for two scalars or if scalar 
%          b = \pm 1 at the moment
function result = mpower(this, varargin)
    % This is a function which involves a second instance of a similar object,
    % so we check if this second instance was also provided
    if length(varargin) ~= 1
        error('Wrong number of arguments in sgem::mpower');
    end
    
    if (numel(varargin{1}) > 1) || ((numel(this) > 1) && (abs(varargin{1}) ~= 1))
        error('Only powers +1 and -1 are supported for the moment')
    end
    
    if numel(this) == 1
        result = this.^varargin{1};
        return;
    end
    
    if varargin{1} == -1
        result = inv(this);
    else
        result = this;
    end
end

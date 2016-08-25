% objectIdentifier - returns the private object identifier. Only friend
%   classes can access this object.
function result = objectIdentifier(this)
    % The sgem class is allowed to access this private
    [ST I] = dbstack('-completenames');
    if (length(ST) < 2) || (isempty(strfind(ST(2).file,'/@gem/')) && isempty(strfind(ST(2).file,'/@sgem/')))
        error('Only gem.m and sgem.m are allowed to access this property.');
    end
    result = this.objectIdentifier;
end

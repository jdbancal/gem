% checkIdentifierValidity - checks whether the c++ pointer class to the
% underlying c++ object that this class instance refers to exists in memory
function result = checkIdentifierValidity(this)
    result = sgem_mex('isValid', this.objectIdentifier);
end

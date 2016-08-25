% precision of each real matrix element
function result = precision(this)
    result = sgem_mex('precision', this.objectIdentifier);
end

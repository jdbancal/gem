% precision of each real matrix element
function result = precision(this)
    result = gem_mex('precision', this.objectIdentifier);
end

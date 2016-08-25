% Convert the data to a table of matlab double
function result = double(this)
    result = gem_mex('double', this.objectIdentifier);
end

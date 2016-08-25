% Convert the data to a sparse table of matlab double
function result = double(this)
    result = sgem_mex('double', this.objectIdentifier);
end

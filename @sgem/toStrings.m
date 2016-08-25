% Convert the data to a matlab strings (or a cell thereof)
function result = toStrings(this)
    result = sgem_mex('toStrings', this.objectIdentifier);
    if (size(result,1) == 1) && (numel(this) == 1)
        result = result{3};
    end
end

% Convert the data to a matlab strings (or a cell thereof)
function result = toStrings(this)
    result = gem_mex('toStrings', this.objectIdentifier);
    if numel(result) == 1
        result = result{1};
    end
end

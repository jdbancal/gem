% Convert the data to a matlab strings (or a cell thereof)
function result = toStrings(this, precision)
    
    if nargin < 2
        precision = -1;
    else
        precision = double(precision);
    end
    
    result = sgem_mex('toStrings', this.objectIdentifier, precision);

    if (size(result,1) == 1) && (numel(this) == 1)
        result = result{3};
    end
end

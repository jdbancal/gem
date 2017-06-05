% isnan - True for NaN elements
function result = isnan(this)

result = logical(sgem_mex('isnan', objectIdentifier(this)));

end

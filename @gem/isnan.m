% isnan - True for NaN elements
function result = isnan(this)

result = logical(gem_mex('isnan', objectIdentifier(this)));

end

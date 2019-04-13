% isinf - True for +Inf and -Inf elements
function result = isinf(this)

result = logical(gem_mex('isinf', this.objectIdentifier));

end

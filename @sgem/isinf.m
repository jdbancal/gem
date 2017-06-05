% isinf - True for +Inf and -Inf elements
function result = isinf(this)

result = logical(sgem_mex('isinf', objectIdentifier(this)));

end

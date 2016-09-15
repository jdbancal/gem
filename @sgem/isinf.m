% isinf - True for +Inf and -Inf elements
function result = isinf(this)

result = (this == Inf) + (this == -Inf);

end

% isfinite - True for all elements except NaN, +Inf and -Inf
function result = isfinite(this)

result = 1 - isnan(this) - isinf(this);

end

% isnan - True for NaN elements
function result = isnan(this)

result = 1 - max((this < 10), (this > -10));

end

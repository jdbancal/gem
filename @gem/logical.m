% logical - returns true for nonzero elements
function result = logical(this)

if sum(isnan(this)) > 0
    error('NaN''s cannot be converted to logicals.');
end

result = (this ~= 0);

end

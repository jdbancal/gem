% isempty - Tells whether the matrix contains no element
function result = isempty(this)

result = (prod(size(this))==0);

end
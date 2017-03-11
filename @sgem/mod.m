% mod - modulo: remainder after division
function result = mod(a, b)

result = a - b.*floor(a./b);

end

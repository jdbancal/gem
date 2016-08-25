% end - returns the last index
function result = end(this, K, N)

if N == 1
    % Then the indexing is done with only one index, as in a(end)
    result = prod(size(this));
else
    % Then the indexing is done with two indices, as in a(end,1:2)
    if (K < 1) || (K > 2)
        error('Gem objects only have two dimensions');
    end
    s = size(this);
    result = s(K);
end

end
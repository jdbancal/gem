% uplus - unary plus operator, doesn't do anything
function result = uplus(this)
    % We just return the object to the output (we don't even need
    % to make a copy, this will be done automatically as soon as
    % the new object (or the old one) will be modified)
    result = this;
end

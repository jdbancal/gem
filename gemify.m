% gemify - turns a matlab dense/sparse matrix into a dense/sparse gem 
function result = gemify(this)
    % If the object is already a gem/sgem, then there is nothing to be done
    if isequal(class(this), 'gem') || isequal(class(this), 'sgem')
        result = this;
        return;
    end
    
    % We create the right object depending on whether the argument is
    % sparse or not
    if issparse(this)
        result = sgem(this);
    else
        result = gem(this);
    end
end

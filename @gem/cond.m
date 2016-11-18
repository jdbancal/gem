% cond - the matrix condition number
function result = cond(this)
    % We compute the ratio of maximum and minumum singular values
    svdMax = svds(this,1);
    svdMin = svds(this,1,'smallest');
    result = svdMax/svdMin;
end

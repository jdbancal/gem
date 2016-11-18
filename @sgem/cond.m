% cond - the matrix condition number
function result = cond(this)
    % We compute the ratio of maximum and minumum singular values
    svdMax = svds(this,1);
    svdMin = svds(this,1,'smallest');
    if (svdMin > svdMax)
        error('minimum sigular value is larger than the maximum one (!)');
    end
    result = svdMax/svdMin;
end

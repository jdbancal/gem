% likeMatlab = gemSparseLikeMatlab(likeMatlab)
% 
% Getting and setting whether the library should deal with sparse matrices
% the way matlab does.
%
% Sometimes, operations on sparse matrices can produce non-sparse results.
% For some of these operations, matlab produces a full matrix. For instance,
%   sparse([0 1 0]) + 2
% gives [2 3 2] in full form. It makes sense to produce a full result here
% because the null element '0' is never preserved upon addition by a non-zero
% number.
% 
% On the other hand, some operations always preserve the zeros elements. It
% then makes sense to encode the result of the operation as a sparse matrix.
% Hence, sparse([0 1 0]) * 2, for instance, gives as a result sparse([0 2 0]).
%
% However, matlab does not always encode results expected to be full in full 
% matrices. For instance, exp(sparse([0 1 0])) returns by default a sparse 
% result, even though 0 is not preserved by the exponentiation : exp(0) = 1.
% This typically results in a full matrix being encoded as a sparse one 
% (which is inefficient).
%
% By default, the gem library keeps the sparse attribute of a matrix after
% an operation only if there is a chance that this result is sparse. Since
% both behaviors may be interesting, this function allows one to define 
% whether we want the gem library to produce sparse or full matrices based 
% on whether the applied function preserves the null element, or whether to
% choose matrix types the way matlab does.
%
% To get the current setting, use
%   gemSparseLikeMatlab
% To set the library to behave like Matlab, use
%   gemSparseLikeMatlab(1)
function likeMatlab = gemSparseLikeMatlab(likeMatlab)
    tmp = sgem(0);
    if nargin < 1
        likeMatlab = getSparseLikeMatlab(tmp);
    else
        setSparseLikeMatlab(tmp,likeMatlab);
    end
end

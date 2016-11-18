% eig - eigenvalues and eigenvectors
%
% At the moment, this function is not implemented
function [V D] = eig(this, varargin)
    error('Eig is not implemented for sparse matrices. Use eigs(x) or eig(full(x)) instead.');
end

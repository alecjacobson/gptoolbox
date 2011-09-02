% RES = innerProd(MTX)
%
% Compute (MTX' * MTX) efficiently (i.e., without copying the matrix)
%
% NOTE: This function used to call a MEX function (C code) to avoid copying, but
% newer versions of matlab have eliminated the overhead of the
% simpler form below. 

function res = innerProd(mtx)

res = mtx' * mtx;

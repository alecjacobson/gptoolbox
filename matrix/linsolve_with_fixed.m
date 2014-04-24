function x = linesolve_with_fixed(A,b,known,c)
% LINSOLVE_WITH_FIXED solves linear system Ax=b with
% known coefficients x(known) 
%
% x = linesolve_with_fixed(A,b,known,c)
%
% Inputs:
%   A  n by n matrix of quadratic coefficients
%   b  n by 1 column of linear coefficients
%   known c by 1 column of indices of known coefficients
%   c  c by 1 column of known coefficients
% Outputs:
%   x  n by cols solution

Nk = length(known);
Uk = size(A,1) - Nk;

% create permutation vectors
idx = 1:size(A,1);
idx(known) = 0;
P = [find(idx');known];
[y,PI] = sort(P);

% reordering (unknowns;knowns)
A = A(P,P);
b = b(P);

% split in unknowns and knowns
U = A(1:Uk,1:Uk);
K = A(1:Uk,Uk+1:end);
b = b(1:Uk);

% solve system
x = U\(b-K*c);

% reorder solution
x = [x;c];
x = x(PI);

end
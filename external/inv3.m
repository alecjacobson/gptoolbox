function X = inv3(A, B, varargin)
% function X = inv3(A, B)
%
% Compute in one shot the matrix-inverse or inversion of multiples (3 x 3)
% linear systems
%
% INPUT:
%   A: (3 x 3 x n) array
%   B: (optional) (3 x p x n) array
% OUTPUT:
%   X:
%   if B is not provided X is (3 x 3 x n) array, X(:,:,k) = inv(A(:,:,k));
%   otherwise, X is (3 x p x n) array, X(:,:,k) = A(:,:,k) \ B(:,:,k).
%
% See also: inv2, eig3, mtimesx
%
% Note: Require James Tursa's MTIMESX
%       http://www.mathworks.com/matlabcentral/fileexchange/25977
%       This function is faster but less accurate than backslash
%       operator when the matrix is ill-conditioned 
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 27-May-2010
%     

if size(A,1) ~= 3 || size(A,2) ~= 3
    error('A must be [3x3xn] array');
end

A = reshape(A, 9, []).';

% minors
m1 = (A(:,5).*A(:,9) - A(:,8).*A(:,6));
m2 = (A(:,2).*A(:,9) - A(:,8).*A(:,3));
m3 = (A(:,2).*A(:,6) - A(:,5).*A(:,3));

X = [+m1 ...
     -m2 ...
     +m3 ...
     -(A(:,4).*A(:,9)-A(:,7).*A(:,6)) ...
     +(A(:,1).*A(:,9)-A(:,7).*A(:,3)) ...
     -(A(:,1).*A(:,6)-A(:,4).*A(:,3)) ...
     +(A(:,4).*A(:,8)-A(:,7).*A(:,5)) ...
     -(A(:,1).*A(:,8)-A(:,7).*A(:,2)) ...
     +(A(:,1).*A(:,5)-A(:,4).*A(:,2))];
     
% Determinant
D = A(:,1).*m1 - A(:,4).*m2 + A(:,7).*m3;

X =  bsxfun(@times, X, 1./D);

X = reshape(X.', 3, 3, []);

if nargin>=2
    X = mtimesx(X,B,'speed');
end

end


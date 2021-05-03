function X = inv2(A, B, varargin)
% function X = inv2(A, B)
%
% Compute in one shot the matrix-inverse or inversion of multiples (2 x 2)
% linear systems
%
% INPUT:
%   A: (2 x 2 x n) array
%   B: (optional) (2 x p x n) array
% OUTPUT:
%   X:
%   if B is not provided X is (2 x 2 x n) array, X(:,:,k) = inv(A(:,:,k));
%   otherwise, X is (2 x p x n) array, X(:,:,k) = A(:,:,k) \ B(:,:,k).
%
% See also: inv3, eig3, mtimesx
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

if size(A,1) ~= 2 || size(A,2) ~= 2
    error('A must be [2x2xn] array');
end

A = reshape(A, 4, []).';

X = [A(:,4) -A(:,2) -A(:,3) A(:,1)];
     
% Determinant
D = A(:,1).*A(:,4) - A(:,2).*A(:,3);

X =  bsxfun(@times, X, 1./D);

X = reshape(X.', 2, 2, []);

if nargin>=2
    X = mtimesx(X,B,'speed');
end

end


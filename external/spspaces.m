function [SpLeft, SpRight] = spspaces(A,opt,tol)
%  PURPOSE: finds left and right null and range space of a sparse matrix A
%
% ---------------------------------------------------
%  USAGE: [SpLeft, SpRight] = spspaces(A,opt,tol)
%
%  INPUT: 
%       A                           a sparse matrix
%       opt                         spaces to calculate
%                                   = 1: left null and range space
%                                   = 2: right null and range space
%                                   = 3: both left and right spaces
%       tol                         uses the tolerance tol when calculating
%                                   null subspaces (optional)
%
%   OUTPUT:
%       SpLeft                      1x4 cell. SpLeft = {} if opt =2.
%           SpLeft{1}               an invertible matrix Q
%           SpLeft{2}               indices, I, of rows of the matrix Q that
%                                   span the left range of the matrix A
%           SpLeft{3}               indices, J, of rows of the matrix Q that
%                                   span the left null space of the matrix A
%                                   Q(J,:)A = 0
%           SpLeft{4}               inverse of the matrix Q
%       SpRight                     1x4 cell. SpRight = {} if opt =1.
%           SpLeft{1}               an invertible matrix Q
%           SpLeft{2}               indices, I, of rows of the matrix Q that
%                                   span the right range of the matrix A
%           SpLeft{3}               indices, J, of rows of the matrix Q that
%                                   span the right null space of the matrix A
%                                   AQ(:,J) = 0
%           SpLeft{4}               inverse of the matrix Q
%
%   COMMENTS:
%       uses luq routine, that finds matrices L, U, Q such that
%
%           A = L | U 0 | Q
%                 | 0 0 |
%       
%       where L, Q, U are invertible matrices, U is upper triangular. This
%       decomposition is calculated using lu decomposition.
%
%       This routine is fast, but can deliver inaccurate null and range
%       spaces if zero and nonzero singular values of the matrix A are not
%       well separated.
%
%   WARNING:
%       right null and rang space may be very inaccurate
%
% Copyright  (c) Pawel Kowal (2006)
% All rights reserved
% LREM_SOLVE toolbox is available free for noncommercial academic use only.
% pkowal3@sgh.waw.pl

if nargin<3
    tol                 = max(max(size(A)) * norm(A,1) * eps,100*eps);
end

switch opt
    case 1
        calc_left       = 1;
        calc_right      = 0;
    case 2
        calc_left       = 0;
        calc_right      = 1;
    case 3
        calc_left       = 1;
        calc_right      = 1;
end

[L,U,Q]                 = luq(A,0,tol);

if calc_left
    if ~isempty(L)
        LL              = L^-1;
    else
        LL              = L;
    end
    S                   = max(abs(U),[],2);
    I                   = find(S>tol);
    if ~isempty(S)
        J               = find(S<=tol);
    else
        J               = (1:size(S,1))';
    end    
    SpLeft              = {LL,I,J,L};
else
    SpLeft              = {};
end
if calc_right
    if ~isempty(Q)
        QQ              = Q^-1;
    else
        QQ              = Q;
    end    
    S                   = max(abs(U),[],1);
    I                   = find(S>tol);
    if ~isempty(S)
        J               = find(S<=tol);
    else
        J               = (1:size(S,2))';
    end
    SpRight             = {QQ,I,J,Q};
else
    SpRight             = {};
end
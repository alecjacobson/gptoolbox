function [L,U,Q] = luq(A,do_pivot,tol)
%  PURPOSE: calculates the following decomposition
%             
%       A = L |Ubar  0 | Q
%             |0     0 |
%
%       where Ubar is a square invertible matrix
%       and matrices L, Q are invertible.
%
% ---------------------------------------------------
%  USAGE: [L,U,Q] = luq(A,do_pivot,tol)
%  INPUT: 
%         A             a sparse matrix
%         do_pivot      = 1 with column pivoting
%                       = 0 without column pivoting
%         tol           uses the tolerance tol in separating zero and
%                       nonzero values
%
%   OUTPUT:
%         L,U,Q          matrices
%
%   COMMENTS:
%         based on lu decomposition
%
% Copyright  (c) Pawel Kowal (2006)
% All rights reserved
% LREM_SOLVE toolbox is available free for noncommercial academic use only.
% pkowal3@sgh.waw.pl

[n,m]                   = size(A);

if ~issparse(A)
    A                   = sparse(A);
end

%--------------------------------------------------------------------------
%       SPECIAL CASES
%--------------------------------------------------------------------------
if size(A,1)==0
    L                   = speye(n);
    U                   = A;
    Q                   = speye(m);
    return;
end
if size(A,2)==0
    L                   = speye(n);
    U                   = A;    
    Q                   = speye(m);
    return;
end        

%--------------------------------------------------------------------------
%       LU DECOMPOSITION
%--------------------------------------------------------------------------
if do_pivot
    [L,U,P,Q]           = lu(A);   
    Q                   = Q';
else
    [L,U,P]             = lu(A);   
    Q                   = speye(m);
end
p                       = size(A,1)-size(L,2);
LL                      = [sparse(n-p,p);speye(p)];
L                       = [P'*L P(n-p+1:n,:)'];
U                       = [U;sparse(p,m)];

%--------------------------------------------------------------------------
%       FINDS ROWS WITH ZERO AND NONZERO ELEMENTS ON THE DIAGONAL
%--------------------------------------------------------------------------
if size(U,1)==1 || size(U,2)==1
    S                   = U(1,1);
else
    S                   = diag(U);
end
I                       = find(abs(S)>tol);
Jl                      = (1:n)';
Jl(I)                   = [];
Jq                      = (1:m)';
Jq(I)                   = [];

Ubar1                   = U(I,I);
Ubar2                   = U(Jl,Jq);
Qbar1                   = Q(I,:);
Lbar1                   = L(:,I);

%--------------------------------------------------------------------------
%       ELININATES NONZEZO ELEMENTS BELOW AND ON THE RIGHT OF THE
%       INVERTIBLE BLOCK OF THE MATRIX U
%
%       UPDATES MATRICES L, Q
%--------------------------------------------------------------------------
if ~isempty(I)
    Utmp                = U(I,Jq);
    X                   = Ubar1'\U(Jl,I)';
    Ubar2               = Ubar2-X'*Utmp;
    Lbar1               = Lbar1+L(:,Jl)*X';

    X                   = Ubar1\Utmp;
    Qbar1               = Qbar1+X*Q(Jq,:);    
    Utmp                = [];
    X                   = [];
end

%Alec:
%LL = [Lbar1 sparse(size(Lbar1,1),size(L,2)-size(Lbar1,2))];
%UU = [Ubar1                                       sparse(size(Ubar1,1),size(U,2)-size(Ubar1,2)); ...
%      sparse(size(U,1)-size(Ubar1,1),size(Ubar1,2)) Ubar2];
%QQ = [Qbar1 ; sparse(size(Q,1)-size(Qbar1,1),size(Q,2)) ];
%norm(full(A-LL*UU*QQ))

%--------------------------------------------------------------------------
%       FINDS ROWS AND COLUMNS WITH ONLY ZERO ELEMENTS
%--------------------------------------------------------------------------
I2                      = find(max(abs(Ubar2),[],2)>tol);
I5                      = find(max(abs(Ubar2),[],1)>tol);

% Alec: These are the cols/rows with zero diagonal but some non-zero elements
I3                      = Jl(I2);
I4                      = Jq(I5);
% Alec: These are the rows/cols (with zero diagonal and) only zero elements
Jq(I5)                  = [];
Jl(I2)                  = [];
U                       = [];

%--------------------------------------------------------------------------
%       FINDS A PART OF THE MATRIX U WHICH IS NOT IN THE REQIRED FORM
%--------------------------------------------------------------------------
% Alec: block corresponding to zero-diagonal but *non-zeros* on off digaonal
A                       = Ubar2(I2,I5);

%--------------------------------------------------------------------------
%       PERFORMS LUQ DECOMPOSITION OF THE MATRIX A
%--------------------------------------------------------------------------
[L1,U1,Q1]              = luq(A,do_pivot,tol);

%--------------------------------------------------------------------------
%       UPDATES MATRICES L, U, Q
%--------------------------------------------------------------------------
Lbar2                   = L(:,I3)*L1;
Qbar2                   = Q1*Q(I4,:);
L                       = [Lbar1 Lbar2 L(:,Jl)];
Q                       = [Qbar1; Qbar2; Q(Jq,:)];

n1                      = length(I);
n2                      = length(I3);
m2                      = length(I4);
U                       = [Ubar1 sparse(n1,m-n1);sparse(n2,n1) U1 sparse(n2,m-n1-m2);sparse(n-n1-n2,m)];
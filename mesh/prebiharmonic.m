function D = prebiharmonic( W,L,a,t,k)
% PREBIHARMONIC Compute the pre-biharmonic kernel SGP 2011 - Raif Rustamov -
% Multiscale Biharmonic Kernels
%
% D = prebiharmonic(W,L,t,k)
%
% Input:
%   W  #V by #V sparse Area weight per vertex (see massmatrix)
%   L  #V by #V sparse Cotan matrix (see cotmatrix)
%   t  Size of support (0..1)
%   a  Total area of the surface
%   k  source vertex
% Outputs:
%   D  #V list of kernel values
%
% Example:
%   [V,F] = readOBJ('woody.obj');
%   W = massmatrix(V,F,'voronoi');
%   L = cotmatrix(V,F);
%   a = sum(doublearea(V,F))/2.0;
%   d = prebiharmonic(W,L,a,0.5,1);



% Scale the t parameter
t = a*t;

% Number of vertices
n = size(W,1);

% x = [f,r]

%% Energy
% 1/2 x'Qx + c'x

% Q = | 0    0  |
%     | 0 2W^-1 |

[i,j,v] = find(W);
i = i+n;
j = j+n;
v = 2.*(ones(n,1)./v);
Q = sparse(i,j,v,2*n,2*n);

c = zeros(2*n,1);

%% Constraints
% lc <= Ax <= uc
% lx <= x  <= ux

% Init
%[~,~,w] = find(2.*ones(n,1)./v);
[~,~,w] = find(W);
Ai = []; Aj = []; Av = [];
lc = []; uc = [];
cur = 1;
lx = -Inf(2*n,1); ux = Inf(2*n,1);
z0 = zeros(n,1); 

% (a) w'f <= t --- single row
[i,j,v] = find([w' z0']);
Ai = [Ai; i'+cur-1];
Aj = [Aj; j'];
Av = [Av; v'];
lc = [lc; -inf];
uc = [uc; t];
cur = cur+1;

% (b) lx(k)=ux(h)=1
lx(k) = 1;
ux(k) = 1;

% (c) Lf-r=0 --- n rows
% becomes [L -eye(n,n)]

%L
[i,j,v] = find(L);
Ai = [Ai; i+cur-1];
Aj = [Aj; j];
Av = [Av; v];

% -eye(n,n)
i = (1:n)';
j = ((1:n)+n)';
v = -ones(n,1);
Ai = [Ai; i+cur-1];
Aj = [Aj; j];
Av = [Av; v];

lc = [lc; zeros(n,1)];
uc = [uc; zeros(n,1)];
cur = cur+n;

% (d) f >= 0 --- n rows
%[i,j,v] = find([eye(n,n), zeros(n,n)]);
i = (1:n)';
j = (1:n)';
v = ones(n,1);

Ai = [Ai; i+cur-1];
Aj = [Aj; j];
Av = [Av; v];
lc = [lc; zeros(n,1)];
uc = [uc; Inf(n,1)];
cur = cur+n;

% build the sparse matrix A
A = sparse(Ai,Aj,Av,cur-1,n*2);

%% Optimize the QP problem
param.MSK_IPAR_LOG = 0;

[res] = mskqpopt(Q,c,A,lc,uc,lx,ux,param);
D = res.sol.itr.xx(1:n);


end


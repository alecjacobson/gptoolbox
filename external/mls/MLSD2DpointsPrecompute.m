function mlsd = MLSD2DpointsPrecompute(varargin)
% MLSD2DPOINTSPRECOMPUTE  Generates a 2D MLSD structure.
%
% function mlsd = MLSD2DpointsPrecompute(p,v,type,a)
%
%  This function precomputes a set of values and generate an MLSD structure
% that embeds infos about the deformation algorithm to be used and embeds
% the parameters that must be used for the process. This function allows to
% initialize a deformation process of one or more an images. For more
% informations see [1].
%
%  [1] "Image Deformation Using Moving Least Squares",
%      Scott Schaefer, Travis McPhail, Joe Warren
%
%  Parameters
%  ----------
% IN:
%  p    = The starting position of the handles points.
%  v    = The points to be moved (i.e. a grid).
%  type = The algorithm type ('affine','similar','rigid'). (def='rigid')
%  a    = The order. (def=2)
%  w    = provide your own weights

% Parsing parameters:
[p,v,type,a,w] = ParseParams(varargin{:});

% Precomputing the weights:
if(isempty(w))
  w = PrecomputeWeights(p,v,a);
end

% Preparing the structure:
mlsd.p = p;
mlsd.v = v;
mlsd.type = type;
mlsd.w = w;
mlsd.constr = 'points';

% Selecting the generation function:
switch type
    case 'affine'
        mlsd.data = PrecomputeAffine(p,v,w);
    case 'similar'
        mlsd.data = PrecomputeSimilar(p,v,w);
    case 'rigid'
        mlsd.data = PrecomputeRigid(p,v,w);
end

% ------------------------ LOCAL FUNCTIONS ------------------------

% Precomputing the weights:
function w = PrecomputeWeights(p,v,a)

% Preparing the output buffer:
w = zeros(size(p,2),size(v,2));

% Iterating on handles:
for i=1:size(p,2)
    % Precomputing the norms^2:
    norms_2 = sum((repmat(p(:,i),[1,size(v,2)])-v).^2,1);
    
    % Precomputing the weights:
    w(i,:) = 1./(norms_2.^a+1e-8);
end
w = w./repmat(sum(w),size(p,2),1);

% -----------------------------------------------------------------

% Precomputing weighted centroids:
function Pstar = PrecomputeWCentroids(p,w)

% Computing the centroids:
Pstar = (p*w)./repmat(sum(w,1),[size(p,1),1]);

% -----------------------------------------------------------------

% Precomputing the affine deformation:
function data = PrecomputeAffine(p,v,w)

% Computing the centroids:
Pstar = PrecomputeWCentroids(p,w);

% Precomputing the first matrix:
M1 = v - Pstar;

% Allocating the central matrix:
a = zeros(1,size(Pstar,2));
b = a;
d = a;

% Iterating on points:
Phat = cell(1,size(p,2));
for i=1:size(p,2)
    % Computing the hat points:
    Phat{i} = repmat(p(:,i),[1,size(Pstar,2)])-Pstar;
    
    % Computing the matrix elements:
    a = a + w(i,:).*Phat{i}(1,:).^2;
    b = b + w(i,:).*Phat{i}(1,:).*Phat{i}(2,:);
    d = d + w(i,:).*Phat{i}(2,:).^2;
end

% Computing the determinant:
det = a.*d - b.^2;

% Computing the inverse:
Ia = d./det;
Ib = -b./det;
Id = a./det;

% Computing the first product elements:
F1 = [sum(M1.*[Ia;Ib],1);sum(M1.*[Ib;Id],1)];

% Computing the A values:
A = zeros(size(p,2),size(Pstar,2));
for j=1:size(p,2)
    % A single element:
    A(j,:) = sum(F1.*Phat{j},1).*w(j,:);
end

% The data structure:
data.A = A;

% -----------------------------------------------------------------

% Precomputing the Asimilar:
function [A,R1] = PrecomputeA(Pstar,Phat,v,w)

% Allocating:
A = cell(1,numel(Phat));

% Fixed part:
R1 = v - Pstar;
R2 = [R1(2,:);-R1(1,:)];

% Iterating on points:
for i=1:numel(Phat)
    % Precomputing the blocks:
    L1 = Phat{i};
    L2 = [L1(2,:);-L1(1,:)];

    % Computing the values:
    A{i}.a = w(i,:).*sum(L1.*R1,1);
    A{i}.b = w(i,:).*sum(L1.*R2,1);
    A{i}.c = w(i,:).*sum(L2.*R1,1);
    A{i}.d = w(i,:).*sum(L2.*R2,1);
end

% -----------------------------------------------------------------

% Precomputing the similar deformation:
function data = PrecomputeSimilar(p,v,w)

% Computing the centroids:
Pstar = PrecomputeWCentroids(p,w);

% Iterating on points:
Phat = cell(1,size(p,2));
mu = zeros(1,size(Pstar,2));
for i=1:size(p,2)
    % Computing the hat points:
    Phat{i} = repmat(p(:,i),[1,size(Pstar,2)])-Pstar;
    
    % Updating the values of mu:
    mu = mu + w(i,:).*sum(Phat{i}.^2,1);
end

% Computing the matrix A:
A = PrecomputeA(Pstar,Phat,v,w);

% Premultiplying A/mu:
for i=1:numel(A)
    % Managing a single matrix:
    A{i}.a = A{i}.a./mu;
    A{i}.b = A{i}.b./mu;
    A{i}.c = A{i}.c./mu;
    A{i}.d = A{i}.d./mu;
end

% The data structure:
data.A = A;

% -----------------------------------------------------------------

% Precomputing the similar deformation:
function data = PrecomputeRigid(p,v,w)

% Computing the centroids:
Pstar = PrecomputeWCentroids(p,w);

% Iterating on points:
Phat = cell(1,size(p,2));
for i=1:size(p,2)
    % Computing the hat points:
    Phat{i} = repmat(p(:,i),[1,size(Pstar,2)])-Pstar;
end

% Computing the matrix A and v-Pstar:
[A,v_Pstar] = PrecomputeA(Pstar,Phat,v,w);

% The norm of v-Pstar:
normof_v_Pstar = sqrt(sum(v_Pstar.^2,1));

% The data structure:
data.A = A;
data.normof_v_Pstar = normof_v_Pstar;

% -----------------------------------------------------------------

% Parsing of parameters:
function [p,v,type,a,w] = ParseParams(varargin)

% Number of parameters:
if nargin<2 error('Too few parameters'); end
if nargin>5 error('Too many parameters'); end

% Set up variables:
varnames = {'p','v','type','a','w'};
for ind=1:nargin
    eval([varnames{ind} ' = varargin{ind} ;']);
end

% Default parameters:
if nargin<3 type='rigid'; end
if nargin<4 a=2; end
if nargin<5 w=[]; end

% Checking the points p and v:
p = points2dnormalize(p);
v = points2dnormalize(v);

% Checking the value of a:
a = abs(a);
if a==0
    a = 1e-8;
end

% Check type param:
if ~strin(type,{'affine','similar','rigid'})
    error('Unknown algorithm type!');
end

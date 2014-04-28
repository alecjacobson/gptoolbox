function fv = MLSD2DTransform(varargin)
% MLSD2DTRANSFORM  Generates a 2D MLSD structure.
%
% function fv = MLSD2DTransform(mlsd,q)
%
%  This function transform the v points. For more informations see [1].
%
%  [1] "Image Deformation Using Moving Least Squares",
%      Scott Schaefer, Travis McPhail, Joe Warren
%
%  Parameters
%  ----------
% IN:
%  mlsd     = A valid mlsd structure.
%  q        = A valid set of new handles (2|3xN points or 4xN segments).
% OUT:
%  fv       = The new points obtained moving v.

% Parsing parameters:
[mlsd,q] = ParseParams(varargin{:});

% Selecting the type of constrint objects:
switch mlsd.constr
    % Points constrained transformation:
    case 'points'
        % Selecting the type of transformation:
        switch mlsd.type
            case 'affine'
                fv = PointsTransformAffine(mlsd,q);
            case 'similar'
                fv = PointsTransformSimilar(mlsd,q);
            case 'rigid'
                fv = PointsTransformRigid(mlsd,q);
        end
    % Segments constrained transformation:
    case 'lines'
        % Selecting the type of transformation:
        switch mlsd.type
            case 'affine'
                fv = LinesTransformAffine(mlsd,q);
            case 'similar'
                fv = LinesTransformSimilar(mlsd,q);
            case 'rigid'
                fv = LinesTransformRigid(mlsd,q);
        end
end

% ------------------------ LOCAL FUNCTIONS ------------------------

% Precomputing weighted centroids:
function Qstar = PrecomputeWCentroids(p,w)

% Computing the centroids:
Qstar = (p*w)./repmat(sum(w,1),[size(p,1),1]);

% -----------------------------------------------------------------

% Precomputing weighted centroids:
function Qstar = PrecomputeWCentroidsLines(q,W,denStar)

% Initializing the Pstar array:
Qstar = zeros(2,size(W.d00,2));

% Iterating on segments:
for i=1:size(q,2)
    % Auxiliary:
    d00_d01 = W.d00(i,:) + W.d01(i,:);
    d01_d11 = W.d01(i,:) + W.d11(i,:);

    % Updating the centroids:
    Qstar = Qstar + repmat(q(1:2,i),[1,size(d00_d01,2)]).*[d00_d01;d00_d01] + ...
                    repmat(q(3:4,i),[1,size(d00_d01,2)]).*[d01_d11;d01_d11];
end

% Computing the centroids:
Qstar = Qstar./repmat(denStar,[2,1]);

% -----------------------------------------------------------------

% An affine transformation:
function fv = PointsTransformAffine(mlsd,q)

% Computing the weighted centroids:
Qstar = PrecomputeWCentroids(q,mlsd.w);
fv = Qstar;

% Adding the affine parts (without the computed translation):
for j=1:size(q,2)
    % Computing the hat points:
    Qhat = repmat(q(:,j),[1,size(Qstar,2)])-Qstar;
    
    % Updating:
    fv = fv + Qhat.*repmat(mlsd.data.A(j,:),[size(Qhat,1),1]);
end

% -----------------------------------------------------------------

% A similar transformation:
function fv = PointsTransformSimilar(mlsd,q)

% Computing the weighted centroids:
Qstar = PrecomputeWCentroids(q,mlsd.w);
fv = Qstar;

% Adding the affine parts (without the computed translation):
for i=1:size(q,2)
    % Computing the hat points:
    Qhat = repmat(q(:,i),[1,size(Qstar,2)])-Qstar;
    
    % Updating:
    fv = fv + [sum(Qhat.*[mlsd.data.A{i}.a;mlsd.data.A{i}.c],1); ...
               sum(Qhat.*[mlsd.data.A{i}.b;mlsd.data.A{i}.d],1)];
end

% -----------------------------------------------------------------

% A similar transformation:
function fv = PointsTransformRigid(mlsd,q)

% Computing the weighted centroids:
Qstar = PrecomputeWCentroids(q,mlsd.w);

% Computing the first step:
fv2 = zeros(size(Qstar));
for i=1:size(q,2)
    % Computing the hat points:
    Qhat = repmat(q(:,i),[1,size(Qstar,2)])-Qstar;
    
    % Adding a vector:
    fv2 = fv2 + [sum(Qhat.*[mlsd.data.A{i}.a;mlsd.data.A{i}.c],1); ...
                 sum(Qhat.*[mlsd.data.A{i}.b;mlsd.data.A{i}.d],1)];
end

% The norms of fv2:
normof_fv2 = sqrt(sum(fv2.^2,1));

% The normalization factor:
norm_fact = mlsd.data.normof_v_Pstar./normof_fv2;

% Generating the output:
fv = fv2.*repmat(norm_fact,[size(fv2,1),1]) + Qstar;

% -----------------------------------------------------------------

% An affine transformation:
function fv = LinesTransformAffine(mlsd,q)

% Computing the weighted centroids:
Qstar = PrecomputeWCentroidsLines(q,mlsd.W,mlsd.data.denStar);
fv = Qstar;

% Adding the affine parts (without the computed translation):
for j=1:size(q,2)
    % Computing the hat points:
    chat = repmat(q(1:2,j),[1,size(Qstar,2)])-Qstar;
    dhat = repmat(q(3:4,j),[1,size(Qstar,2)])-Qstar;
    
    % Updating:
    fv = fv + [mlsd.data.A{j}.a;mlsd.data.A{j}.a].*chat + ...
              [mlsd.data.A{j}.b;mlsd.data.A{j}.b].*dhat;
end

% -----------------------------------------------------------------

% An affine transformation:
function fv = LinesTransformSimilar(mlsd,q)

% Computing the weighted centroids:
Qstar = PrecomputeWCentroidsLines(q,mlsd.W,mlsd.data.denStar);
fv = Qstar;

% Adding the affine parts (without the computed translation):
for j=1:size(q,2)
    % Computing the hat points:
    chat = repmat(q(1:2,j),[1,size(Qstar,2)])-Qstar;
    dhat = repmat(q(3:4,j),[1,size(Qstar,2)])-Qstar;
    
    % Updating:
    L = [chat;dhat];
    A = mlsd.data.A{j};
    fv = fv + [sum(L.*[A.a11;A.a21;A.a31;A.a41],1);
               sum(L.*[A.a12;A.a22;A.a32;A.a42],1)];
end

% -----------------------------------------------------------------

% An affine transformation:
function fv = LinesTransformRigid(mlsd,q)

% Computing the weighted centroids:
Qstar = PrecomputeWCentroidsLines(q,mlsd.W,mlsd.data.denStar);

% Adding the affine parts (without the computed translation):
fv2 = zeros(size(Qstar));
for j=1:size(q,2)
    % Computing the hat points:
    chat = repmat(q(1:2,j),[1,size(Qstar,2)])-Qstar;
    dhat = repmat(q(3:4,j),[1,size(Qstar,2)])-Qstar;
    
    % Adding a vector:
    L = [chat;dhat];
    A = mlsd.data.A{j};
    fv2 = fv2 + [sum(L.*[A.a11;A.a21;A.a31;A.a41],1);
                 sum(L.*[A.a12;A.a22;A.a32;A.a42],1)];
end

% The norms of fv2:
normof_fv2 = sqrt(sum(fv2.^2,1));
normof_fv2(find(normof_fv2==0)) = 1e-8;

% The normalization factor:
norm_fact = mlsd.data.normof_v_Pstar./normof_fv2;

% Generating the output:
fv = fv2.*repmat(norm_fact,[size(fv2,1),1]) + Qstar;

% -----------------------------------------------------------------

% Parsing of parameters:
function [mlsd,q] = ParseParams(varargin)

% Number of parameters:
if nargin<2 error('Too few parameters'); end
if nargin>2 error('Too many parameters'); end

% Set up variables:
varnames = {'mlsd','q'};
for ind=1:nargin
    eval([varnames{ind} ' = varargin{ind} ;']);
end

% The structure mlsd is not checked, it must be valid.

% Checking the points q:
if size(q,1)<4
    q = points2dnormalize(q);
end

% Checking the shape of q:
if size(q,2)~=size(mlsd.p,2)
    error('The new handle must be the same number!');
end

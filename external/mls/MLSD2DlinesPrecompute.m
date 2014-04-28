function mlsd = MLSD2DlinesPrecompute(varargin)
% MLSD2DLINESPRECOMPUTE  Generates a 2D MLSD structure.
%
% function mlsd = MLSD2DlinesPrecompute(p,v,type)
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
%  p    =  The starting position of the handles lines (as [a;b] vectors with
%         a and b column points that define the extrema of the segment).
%  v    = The points to be moved (i.e. a grid).
%  type = The algorithm type ('affine','similar','rigid'). (def='rigid')

% Parsing parameters:
[p,v,type] = ParseParams(varargin{:});

% Removing coincidences:
v = RemoveCoincidences(v,p);

% Precomputing the weights:
W = PrecomputeWeights(p,v);

% Precomputing the centroids:
[Pstar,denStar] = PrecomputeWCentroids(p,W);

% Preparing the structure:
mlsd.p = p;
mlsd.v = v;
mlsd.type = type;
mlsd.W = W;
mlsd.constr = 'lines';

% Selecting the generation function:
switch type
    case 'affine'
        mlsd.data = PrecomputeAffine(p,v,W,Pstar);
    case 'similar'
        mlsd.data = PrecomputeSimilar(p,v,W,Pstar);
    case 'rigid'
        mlsd.data = PrecomputeRigid(p,v,W,Pstar);
end

% Other data:
mlsd.data.denStar = denStar;

% ------------------------ LOCAL FUNCTIONS ------------------------

% Removing coincidences:
function v = RemoveCoincidences(v,p)

% Iterating on points:
for i=1:size(p,2)
    % Coincidences with a:
    ind = find(v(1,:)==p(1,i) & v(2,:)==p(2,i));
    v(:,ind) = v(:,ind) + 1e-8;
    
    % Coincidences with b:
    ind = find(v(1,:)==p(3,i) & v(2,:)==p(4,i));
    v(:,ind) = v(:,ind) + 1e-8;
end

% -----------------------------------------------------------------

% Precomputing the weights:
function W = PrecomputeWeights(p,v)

% Init W:
W.d00 = zeros(size(p,2),size(v,2));
W.d01 = W.d00;
W.d11 = W.d00;

% Computing usefull quantities:
a_b = p(1:2,:)-p(3:4,:);
b_a = -a_b;
for i=1:size(p,2)
    % Computing the a-v  and b-v values:
    a_v = repmat(p(1:2,i),[1,size(v,2)])-v;
    b_v = repmat(p(3:4,i),[1,size(v,2)])-v;
    v_b = -b_v;
    
    % Computing the triplaned values:
    a_v_t = [-a_v(2,:);a_v(1,:)];
    b_v_t = [-b_v(2,:);b_v(1,:)];
    
    % Computing the big Delta and all the other combinations:
    Num1 = b_v(1,:).*b_a(1,i) + b_v(2,:).*b_a(2,i);
    Num2 = a_v(1,:).*a_b(1,i) + a_v(2,:).*a_b(2,i);
    Den1 = b_v_t(1,:).*b_a(1,i) + b_v_t(2,:).*b_a(2,i);
    Delta = a_v_t(1,:).*a_b(1,i) + a_v_t(2,:).*a_b(2,i);
    
    % Computing the angle theta:
    theta = atan2(Den1,Num1)-atan2(Delta,Num2);
    
    % Clearing:
    clear Num1 Num2 Den1 a_v_t b_v_t;
    
    % The beta vectors:
    beta00 = sum(a_v.^2,1);
    beta01 = sum(a_v.*v_b,1);
    beta11 = sum(v_b.^2,1);
    
    %% Precomputing values:
    %a_bNorm = norm(a_b(:,i));
    %a_bNorm_2Delta2 = a_bNorm./(2*Delta.^2);
    %theta_Delta = theta./Delta;
    %
    %% Computing the weights:
    %W.d00(i,:) = a_bNorm_2Delta2.*(beta01./beta00-beta11.*theta_Delta);
    %W.d01(i,:) = a_bNorm_2Delta2.*(1-beta01.*theta_Delta);
    %W.d11(i,:) = a_bNorm_2Delta2.*(beta01./beta11-beta00.*theta_Delta);

    % Getting the case vector:
    DeltaIsZero = Delta==0;
    
    % The non-zero case:
    ind = find(~DeltaIsZero);
    a_bNorm = norm(a_b(:,i));
    a_bNorm_2Delta2 = a_bNorm./(2*Delta(ind).^2);
    theta_Delta = theta(ind)./Delta(ind);
    W.d00(i,ind) = a_bNorm_2Delta2.*(beta01(ind)./beta00(ind)-beta11(ind).*theta_Delta);
    W.d01(i,ind) = a_bNorm_2Delta2.*(1-beta01(ind).*theta_Delta);
    W.d11(i,ind) = a_bNorm_2Delta2.*(beta01(ind)./beta11(ind)-beta00(ind).*theta_Delta);
    
    
    % The zero case:
    ind = find(DeltaIsZero);
    if ~isempty(ind)
        % Computing the other elements:
        a_bNorm5 = a_bNorm.^5;
        a_vb_a = a_v(1,ind).*b_a(1,i) + a_v(2,ind).*b_a(2,i);
        v_bb_a = v_b(1,ind).*b_a(1,i) + v_b(2,ind).*b_a(2,i);
        W.d00(i,ind) = a_bNorm5./(3*v_bb_a.*a_vb_a.^3+1e-8);
        W.d01(i,ind) = a_bNorm5./(6*v_bb_a.^2.*a_vb_a.^2+1e-8);
        W.d11(i,ind) = a_bNorm5./(3*v_bb_a.^3.*a_vb_a+1e-8);
    end
    
    % Clearing:
    clear a_vb_a v_bb_a ind DeltaIsZero a_v b_v v_b;
    % Clearing:
    clear beta00 beta01 beta11 theta Delta theta_Delta;
end

% -----------------------------------------------------------------

% Precomputing weighted centroids:
function [Pstar,denStar] = PrecomputeWCentroids(p,W)

% Computing the den star:
denStar = sum(W.d00+2*W.d01+W.d11,1);

% Initializing the Pstar array:
Pstar = zeros(2,size(W.d00,2));

% Iterating on segments:
for i=1:size(p,2)
    % Auxiliary:
    d00_d01 = W.d00(i,:) + W.d01(i,:);
    d01_d11 = W.d01(i,:) + W.d11(i,:);

    % Updating the centroids:
    Pstar = Pstar + repmat(p(1:2,i),[1,size(d00_d01,2)]).*[d00_d01;d00_d01] + ...
                    repmat(p(3:4,i),[1,size(d00_d01,2)]).*[d01_d11;d01_d11];
end

% Computing the centroids:
Pstar = Pstar./[denStar;denStar];

% -----------------------------------------------------------------

% Precomputing the affine deformation:
function data = PrecomputeAffine(p,v,W,Pstar)

% Computing the ahat and bhat and the matrix elements:
ahat = cell(1,size(p,2));
bhat = ahat;
a = zeros(1,size(Pstar,2));
b = a;
d = a;
for i=1:size(p,2)
    % Computing ahat{i} and  bhat{i}:
    ahat{i} = repmat(p(1:2,i),[1,size(Pstar,2)])-Pstar;
    bhat{i} = repmat(p(3:4,i),[1,size(Pstar,2)])-Pstar;
    
    % Computing the matrix elements:
    a = a + ahat{i}(1,:).^2.*W.d00(i,:) + ...
            2*ahat{i}(1,:).*bhat{i}(1,:).*W.d01(i,:) + ...
            bhat{i}(1,:).^2.*W.d11(i,:);
	b = b + ahat{i}(1,:).*ahat{i}(2,:).*W.d00(i,:) + ...
            (ahat{i}(2,:).*bhat{i}(1,:)+ahat{i}(1,:).*bhat{i}(2,:)).*W.d01(i,:) + ...
            bhat{i}(1,:).*bhat{i}(2,:).*W.d11(i,:);
    d = d + ahat{i}(2,:).^2.*W.d00(i,:) + ...
            2*ahat{i}(2,:).*bhat{i}(2,:).*W.d01(i,:) + ...
            bhat{i}(2,:).^2.*W.d11(i,:);
end

% Computing the inverse:
det = a.*d - b.^2;
det(find(det==0)) = 1e-8;
Ia = d./det;
Ib = -b./det;
Id = a./det;

% Computing the v-Pstar values:
M1 = v - Pstar;

% Computing the first product elements:
F1 = [sum(M1.*[Ia;Ib],1);sum(M1.*[Ib;Id],1)];

% Computing the A values:
A = cell(1,size(p,2));
for j=1:size(p,2)
    % All but W:
    a = sum(F1.*ahat{j},1);
    b = sum(F1.*bhat{j},1);
    
    % Multipling per W:
    A{j}.a = a.*W.d00(j,:) + b.*W.d01(j,:);
    A{j}.b = a.*W.d01(j,:) + b.*W.d11(j,:);
end

% The data structure:
data.A = A;

% -----------------------------------------------------------------

% Precomputing the Asimilar:
function [A,R1] = PrecomputeA(Pstar,ahat,bhat,v,W)

% Allocating:
A = cell(1,numel(ahat));

% Fixed part:
R1 = v - Pstar;
R2 = [R1(2,:);-R1(1,:)];
% Iterating on segments:
for i=1:numel(ahat)
    % The transformed versions:
    ahat_t = [ahat{i}(2,:);-ahat{i}(1,:)];
    bhat_t = [bhat{i}(2,:);-bhat{i}(1,:)];

    % Computing the first block:
    a11 = sum(ahat{i}.*R1,1);
    a12 = sum(ahat{i}.*R2,1);
    a21 = sum(ahat_t.*R1,1);
    a22 = sum(ahat_t.*R2,1);
    a31 = sum(bhat{i}.*R1,1);
    a32 = sum(bhat{i}.*R2,1);
    a41 = sum(bhat_t.*R1,1);
    a42 = sum(bhat_t.*R2,1);
    
    % Computing the values:
    A{i}.a11 = a11.*W.d00(i,:) + a31.*W.d01(i,:);
    A{i}.a12 = a12.*W.d00(i,:) + a32.*W.d01(i,:);
    A{i}.a21 = a21.*W.d00(i,:) + a41.*W.d01(i,:);
    A{i}.a22 = a22.*W.d00(i,:) + a42.*W.d01(i,:);
    A{i}.a31 = a11.*W.d01(i,:) + a31.*W.d11(i,:);
    A{i}.a32 = a12.*W.d01(i,:) + a32.*W.d11(i,:);
    A{i}.a41 = a21.*W.d01(i,:) + a41.*W.d11(i,:);
    A{i}.a42 = a22.*W.d01(i,:) + a42.*W.d11(i,:);
end

% -----------------------------------------------------------------

% Precomputing the similar deformation:
function data = PrecomputeSimilar(p,v,W,Pstar)

% Iterating on points:
ahat = cell(1,size(p,2));
bhat = ahat;
mu = zeros(1,size(v,2));
for i=1:size(p,2)
    % Computing ahat{i} and  bhat{i}:
    ahat{i} = repmat(p(1:2,i),[1,size(Pstar,2)])-Pstar;
    bhat{i} = repmat(p(3:4,i),[1,size(Pstar,2)])-Pstar;
    
    % Computing the mu values:
    mu = mu + sum(ahat{i}.^2,1).*W.d00(i,:) + ...
              2*sum(ahat{i}.*bhat{i},1).*W.d01(i,:) + ...
              sum(bhat{i}.^2,1).*W.d11(i,:);
end

% Computing the matrix A:
A = PrecomputeA(Pstar,ahat,bhat,v,W);

% Premultiplying A/mu:
for i=1:numel(A)
    % Managing a single matrix:
    A{i}.a11 = A{i}.a11./mu;
    A{i}.a12 = A{i}.a12./mu;
    A{i}.a21 = A{i}.a21./mu;
    A{i}.a22 = A{i}.a22./mu;
    A{i}.a31 = A{i}.a31./mu;
    A{i}.a32 = A{i}.a32./mu;
    A{i}.a41 = A{i}.a41./mu;
    A{i}.a42 = A{i}.a42./mu;
end

% The data structure:
data.A = A;

% -----------------------------------------------------------------

% Precomputing the similar deformation:
function data = PrecomputeRigid(p,v,W,Pstar)

% Iterating on points:
ahat = cell(1,size(p,2));
bhat = ahat;
for i=1:size(p,2)
    % Computing ahat{i} and  bhat{i}:
    ahat{i} = repmat(p(1:2,i),[1,size(Pstar,2)])-Pstar;
    bhat{i} = repmat(p(3:4,i),[1,size(Pstar,2)])-Pstar;
end

% Computing the matrix A and v-Pstar:
[A,v_Pstar] = PrecomputeA(Pstar,ahat,bhat,v,W);

% The norm of v-Pstar:
normof_v_Pstar = sqrt(sum(v_Pstar.^2,1));

% The data structure:
data.A = A;
data.normof_v_Pstar = normof_v_Pstar;

% -----------------------------------------------------------------

% Parsing of parameters:
function [p,v,type] = ParseParams(varargin)

% Number of parameters:
if nargin<2 error('Too few parameters'); end
if nargin>3 error('Too many parameters'); end

% Set up variables:
varnames = {'p','v','type'};
for ind=1:nargin
    eval([varnames{ind} ' = varargin{ind} ;']);
end

% Default parameters:
if nargin<3 type='rigid'; end

% Checking the points p and v:
v = points2dnormalize(v);

% Check type param:
if ~strin(type,{'affine','similar','rigid'})
    error('Unknown algorithm type!');
end

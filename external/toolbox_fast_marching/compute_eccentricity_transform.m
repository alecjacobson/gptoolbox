function [E,Ind] = compute_eccentricity_transform(M, options)

% compute_eccentricity_transform - compute the eccentricity of a binary shape
%
%   [E,Ind] = compute_eccentricity_transform(M, options);
%
%   M is a binary shape, should be 0 outside the shape, and 1 inside.
%
%   E is the s-eccentriciy of the shape, ie
%       E(i) = (1/p sum_j d(i,j)^s)^(1/s)       if s<0
%       E(i) = max_j d(i,j)                     if s=Inf
%   where d(i,j) is the distance between points and j spans
%   the set of p points on the boundary of M.
%   The value of s is set in options.s.
%
%   If options.metric='euclidean' then d(x,y) is the euclidean distance |x-y|.
%   If options.metric='euclidean' then d(x,y) is the geodesic distance
%       inside the shape.
%
%   To reduce the number of samples used, you can set options.nb_samples (e.g. 300).
%
%   For s==Inf only:
%       Ind(i) is the index of the point j that reaches E(i)=d(i,Ind(i)).
%
%   Copyright (c) 2006 Gabriel Peyr?

options.null = 0;

if size(M,3)>1
    E = [];
    for i=1:size(M,3)
        [E(:,:,i),Ind(:,:,i)] = compute_eccentricity_transform(M(:,:,i), options);
    end
    return;
end

if isfield(options, 'nb_samples')
    nb_samples = options.nb_samples;
else
    nb_samples = 50;
end
if isfield(options, 's')
    s = options.s;
else
    s = Inf;
end
if isfield(options, 'metric')
    metric = options.metric;
else
    metric = 'geodesic';
end

n = size(M,1);
epsilon = 1e-10;


options.nb_iter_max = Inf;
options.Tmax = 2*n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do an initial central propagation in order to have only 1 connected component
W = epsilon+M;
I = find(M(:)==1); x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x); 
[tmp,c] = min( X(I).^2 + Y(I).^2 );
[a b] = ind2sub([n n],I(c));
options.constraint_map = [];
[d,tmp] = perform_fast_marching(W, [a;b], options);
I = find(d>1e5 | d<0); M(I) = 0;

L = zeros(n)-Inf; L(M==1) = +Inf;
options.constraint_map = L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recompute potential
W = ones(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute sampling locations
if s==Inf
    % points on the boundary
    h = [0 1 0; 1 1 1; 0 1 0];
    h = h./sum(h(:));
    Mh = perform_convolution(M,h);
    samples = find(Mh>0 & M==0);
else
    % points inside
    samples = find(M==1);
end


% perform sub-sampling
p = length(samples);
nb_samples = min(nb_samples,p);
sel = randperm(p); sel = sel(1:nb_samples);
samples = samples(sel);
[a b] = ind2sub([n n],samples);

% perform FM for every point on the boudary
E = zeros(n);
Ind = zeros(n,n);
for i=1:nb_samples
    progressbar(i,nb_samples);
    start_points = [a(i);b(i)];
    if strcmp(metric, 'geodesic')
        [d,tmp] = perform_fast_marching(W, start_points, options);
    else
        % compute euclidean distance
        d = compute_euclidean_distance(M, start_points);
    end
    d(M==0) = 0;
    d(d==Inf) = 0;
    d(d>1e5) = 0;
    if s==Inf
        E = max(E,d);
        Ind(E==d) = i;  
    else
        E = E + d.^s;
    end  
end
if s~=Inf
    E = (E/nb_samples).^(1/s);
end

if s==Inf
    Ind = samples(Ind);
end
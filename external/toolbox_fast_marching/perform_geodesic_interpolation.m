function [M,G] = perform_geodesic_interpolation(W,points,f,options)

% perform_geodesic_interpolation - interpolate function values
%
%   [M,G] = perform_geodesic_interpolation(W,points,f,options);
%
%   f(:,i) is the value of the function at the point points(:,i).
%
%   options.method can be 'powerlaw' or 'gaussian'.
%
%   Copyright (c) 2007 Gabriel Peyré


n = size(W,1);
s = size(f,1); % dimenstionality
npoints = size(points,2);
points = round(points);
points = clamp(points, 1,n);

options.null = 1;
if isfield(options,'method')
    method = options.method;
else
    method = 'powerlaw';
end
if isfield(options,'sigma')
    sigma = options.sigma;
else
    sigma = 5/n;
end
if isfield(options,'sigma')
    alpha = options.alpha;
else
    alpha = 3;
end

% compute geodesic distances
G = zeros(n,n,1,npoints);

for i=1:npoints
    progressbar(i,npoints);
    [G(:,:,1,i),Z,Q] = perform_fast_marching(W, points(:,i), options);
end

% perform the interpolation
switch lower(method)
    case 'powerlaw'
        S = 1 ./ ( G.^alpha + sigma^alpha );
    case 'gaussian'
        S = exp( -(G/sigma).^2 );
    otherwise
        error('Unknown method.');
end

D = repmat( sum(S,4), [1 1 1 npoints] );
D(D<eps) = 1;
S = S ./ D;
S = repmat(S, [1 1 s 1]);
T = repmat( reshape(f, [1 1 s npoints]), [n n 1 1] );
M = S .* T;
M = sum( M, 4 );
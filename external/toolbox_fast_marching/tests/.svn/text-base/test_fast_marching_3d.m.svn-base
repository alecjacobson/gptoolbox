% test for fast marching
%
%   Copyright (c) 2004 Gabriel Peyré

n = 40;

% gaussian weight (path will avoid center of the cube)
x = -1:2/(n-1):1;
[X,Y,Z] = meshgrid(x,x,x);
sigma = 0.4;
W = 1./(1 + exp( -(X.^2+Y.^2+Z.^2)/sigma^2 ) );

k = 5;
start_points = [n-k;k;k];
end_points = [k;n-k;n-k];

options.nb_iter_max = Inf;

[D,S] = perform_fast_marching(W, start_points, options);
path = compute_geodesic(D,end_points);

% draw the path
plot_fast_marching_3d(D,S,path,start_points,end_points);
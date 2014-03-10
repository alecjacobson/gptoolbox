% test for circular shape priors
%
%   Copyright (c) 2004 Gabriel Peyré

% some parameters
n = 400;
name = 'diatom';
name = 'surbiscbox';
name = 'stephanodiscusniagarae';

% load test imges
file_name = ['data/', name, '.png'];
M_original = double( imread(file_name) );
M_original = rescale(M_original);
M = M_original;

% compute the speed function
W = compute_edge_energy(M);

% resize
m = n/max(size(M));
m = round( size(M)*m );
warning off;
M = perform_image_resize(M,m(1),m(2));
W = perform_image_resize(W,m(1),m(2));
warning on;

% pick start/end/center points
[start_points,end_points,center_point] = pick_start_end_point(M, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale the metric to avoid shrinkage near zeros
if 0
[Y,X] = meshgrid(1:m(2), 1:m(1));
d = sqrt( (X-center_point(1)).^2 + (Y-center_point(2)).^2 );
d(center_point(1),center_point(2)) = 1e-9;
W = W.*d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute heuristic mapping
[Y,X] = meshgrid(1:m(1), 1:m(2));
r = sqrt( (center_point(1)-start_points(1)).^2 + (center_point(2)-start_points(2)).^2 );
s = sqrt( (X-start_points(1)).^2 + (Y-start_points(2)).^2 );
I = find(s>2*r);
H = r*( 2*pi - 2*asin(s/(2*r)) );
H(I) = pi*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform front propagation
clear options;
options.nb_iter_max = Inf;
options.end_points = end_points;
disp('Performing front propagation.');
[D,S] = perform_circular_fast_marching_2d(W, start_points, center_point, options, H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract the path
disp('Extracting path.');
% end_points = start_points; end_points(2) = end_points(2)-4;
options.Tmax = sum(m)*2;
options.start_points = start_points;
path = compute_geodesic(D,end_points, options);

% litle cleaning in the path ...
if 0
d = (path(:,1)-start_points(1)).^2 + (path(:,2)-start_points(2)).^2;
p = length(d);
[tmp,I] = min(d(round(p/2):end));
path = path(1:(p/2+I(1)), :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display
% scale the path so that it fits into the original image
path = (path-1)*size(M,1)/m(1)+1;
start_points = (start_points-1)*size(M,1)/m(1)+1;
end_points = (end_points-1)*size(M,1)/m(1)+1;
plot_fast_marching_2d(M,[],path,start_points,end_points);
colormap gray(256);
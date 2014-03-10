% test for segmentation using front propagation
%
%   Copyright (c) 2006 Gabriel Peyré

n = 200;
map_name = 'data/room.png';
map_name = 'data/cavern.png';


[W, cm] = imread(map_name);
W = rescale( double(W) );
W = W.^0.2; % enhance contrast
W = W + 1e-6; % put low speed on background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick starting point
start_points = pick_start_end_point(W);

options.nb_iter_max = Inf;
options.Tmax = Inf;
disp('Performing front propagation.');
[D,S] = perform_fast_marching(W, start_points, options);

T = sum(size(W)); % threshold value
clf;
subplot(1,2,1);
imagesc(W);
title('Map');
axis image; axis off;
subplot(1,2,2);
imagesc(D>T);
title('Segmented');
axis image; axis off;
colormap gray(256);


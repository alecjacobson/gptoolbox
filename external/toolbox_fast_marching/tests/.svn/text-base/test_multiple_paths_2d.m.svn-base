% test for multiple path extraction
%
%   Copyright (c) 2005 Gabriel Peyré

n = 256;

% gaussian bump
sigma = 0.5;
x = linspace(-1,1,n);
[Y,X] = meshgrid(x);
W = exp(-(X.^2+Y.^2)/sigma^2);
W = rescale( double(W) );
W = W + 2;   
W = 1./W;

% startpoints/endpoints
start_points = [1:n; ones(1,n)];
N = 20;
end_points = 1 + round( linspace(0,1,N)*(n-1) );
end_points = [end_points; ones(1,N)*n];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using line to point distance function
% perform marching
disp('Performing front propagation.');
clear options;
options.nb_iter_max = Inf;
[D,S] = perform_fast_marching(W, start_points, options);

disp('Extracting path.');
paths1 = compute_geodesic(D,end_points, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using point to point distance function
% perform marching
start_points = [end_points(1,:); ones(1,N)];
paths2 = {};
for i=1:size(end_points,2)
    disp('Performing front propagation.');
    clear options;
    options.nb_iter_max = Inf;
    options.end_points = end_points(:,i);
    [D,S] = perform_fast_marching(W, start_points(:,i), options);
    paths2{end+1} = compute_geodesic(D,end_points(:,i), options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display
clf;
subplot(1,2,1);
plot_fast_marching_2d(W,[],paths1,[],[]);
title('Point-to-line geodesics');
colormap gray(256);

subplot(1,2,2);
plot_fast_marching_2d(W,[],paths2,[],[]);
title('Point-to-point geodesics');
colormap gray(256);

if 0
    clf;
    plot_fast_marching_2d(W,[],paths1,[],[]);
    title('Point-to-line geodesics');
    colormap gray(256);
    saveas(gcf, 'tmp1.eps', 'eps');

    clf;
    plot_fast_marching_2d(W,[],paths2,[],[]);
    title('Point-to-point geodesics');
    colormap gray(256);
    saveas(gcf, 'tmp2.eps', 'eps');
end
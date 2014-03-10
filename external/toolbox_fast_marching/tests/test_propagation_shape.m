% test for propagation and geodesic extraction on 2D planar shape

path(path, 'toolbox/');
path(path, 'data/');

name = 'chicken';
name = 'apple';
name = 'cavern';
name = 'camel';
name = 'giraffe';


rep = 'results/shape-geodesics/';
if not(exist(rep))
    mkdir(rep);
end

n = 128;
M = rescale( load_image(name,n), 0,1 );
M = perform_blurring(M,5);
M = double(M>0.5);

% make sure pixels on the boundary are black
if M(1)==1
    M = 1-M;
end

warning off;
imwrite(1-M, [rep name '-shape.png'], 'png');
warning off;

% compute geodesic distance
clf;
imagesc(M); axis image; axis off;
title('click on a point inside the shape');
[y,x] = ginput(1);
start_points = round([x y]');
W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;
options.constraint_map = L;
disp('Compute distance function');
[D,S,Q] = perform_fast_marching(W, start_points, options);


bound = compute_shape_boundary(M);
nbound = size(bound,1);
npaths = 30;
sel = round(linspace(1,nbound+1,npaths+1)); sel(end) = [];
end_points = bound(sel,:);

disp('Extract paths');
paths = {};
D1 = D; D1(M==0) = 1e9;
for i=1:npaths
    paths{i} = compute_geodesic(D1,end_points(i,:)');
%    paths{i} = compute_discrete_geodesic(D1,end_points(i,:)')';
end

ms = 30; lw = 3;
% display
A = convert_distance_color(D);
clf; hold on;
imageplot(A); axis image; axis off;
for i=1:npaths
    end_point = end_points(i,:);
    h = plot( paths{i}(2,:), paths{i}(1,:), 'k' );
    set(h, 'LineWidth', lw);    
    h = plot(end_point(2),end_point(1), '.b');
    set(h, 'MarkerSize', ms);    
end
h = plot(start_points(2),start_points(1), '.r');
set(h, 'MarkerSize', ms);
hold off;
colormap jet(256);
axis ij;
saveas(gcf, [rep name '-geodesics.png'], 'png');

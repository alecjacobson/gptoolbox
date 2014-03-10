% test for propagation constrained to a sub-set of points


n = 300;
name = 'mountain';
name = 'constant';

% number of points
p = 80;
start_points = floor(rand(2,p)*(n-1))+1;

x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
L = zeros(n)-Inf;  
I = find(X.^2 + Y.^2 < 0.6^2); L(I) = Inf;
options.constraint_map = L;

[M,W] = load_potential_map(name, n);

[D,Z,Q] = perform_fast_marching(W, start_points, options);

clf;
subplot(1,2,1);
hold on;
imagesc(Q'); axis image; axis off;
plot(start_points(1,:), start_points(2,:), 'k.');
hold off;
subplot(1,2,2);
hold on;
imagesc(D'); axis image; axis off;
plot(start_points(1,:), start_points(2,:), 'k.');
hold off;
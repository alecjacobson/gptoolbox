% test for Voronoi cells extractions


n = 300;
name = 'mountain';
name = 'constant';

% number of points
p = 80;
start_points = floor(rand(2,p)*(n-1))+1;

[M,W] = load_potential_map(name, n);

[D,Z,Q] = perform_fast_marching(W, start_points);

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

faces = compute_voronoi_triangulation(Q, vertex);
edges = compute_edges(faces);

hold on;
imagesc(Q'); axis image; axis off;
plot_edges(edges, start_points, 'k');
hold off;
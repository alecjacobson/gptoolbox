% test for meshing the inside of an object

path(path,'toolbox/');
path(path,'data/');
name = 'mm';
n = 400;
Ma = load_image(name,n-10);
Ma = rescale(abs(Ma));

% to avoid that the shape touches the boundary
M = zeros(n,n,3);
M(6:n-5,6:n-5,:) = Ma;

repimg = 'results/meshing-shape/';
if ~exist(repimg)
    mkdir(repimg);
end


M1 = mean(M,3);
mask = 1-(M1==M1(1));


boundary = compute_shape_boundary(mask)';
Ibound = boundary(1,:) + (boundary(2,:)-1)*n;
for k=1:size(M,3)
    Ma = M(:,:,k); 
    Ma(M1<=0) = 1; 
    Ma(Ibound) = 1;
    M(:,:,k) = Ma;
end

% number of samples for the mesh
if not(exist('p'))
    p = 400;
end
% use an adaptive distance field
if not(exist('use_adaptive'))
    use_adaptive = 1;
end

if use_adaptive
    % compute boundary points
    boundary = compute_shape_boundary(mask)';
    % compute distance to boundary
    [D,Z,Q] = perform_fast_marching(ones(n), boundary);
    % set up a distancing field
    R = 0.8;
    D1 = min(rescale(D),R);
    H = sqrt( R^2 - (D1-R).^2 ) * n;
    W = rescale( D, 0.1,1 );
else
    W = ones(n);
end



%% perform sampling using farthest point
L = zeros(n) - Inf;
I = find(mask); L(I) = Inf;
vertex = [n/2;n/2];
options.constraint_map = L;
vertex = perform_farthest_point_sampling(W, vertex, p, options );

%% compute the associated triangulation
[D,Z,Q] = perform_fast_marching(W, vertex, options);
faces = compute_voronoi_triangulation(Q, vertex);

%% display
clf;
hold on;
imagesc(rescale(M)); axis image; axis off;
plot_edges(compute_edges(faces), vertex(2:-1:1,:), 'r');
plot(vertex(2,:), vertex(1,:), 'b.', 'MarkerSize', 8);
hold off;
axis tight; axis image; axis off;
colormap gray(256);
axis ij;

str = [name '-mesh-' num2str(p)];
if use_adaptive
    str = [str '-adaptive'];
end
saveas(gcf, [repimg str '.png'], 'png');
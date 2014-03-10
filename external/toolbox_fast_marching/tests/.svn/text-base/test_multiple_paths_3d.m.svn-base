% test for multiple path extraction in 3D
n = 61;

% create a "vase" 3D function
x = linspace(-1,1,n);
[X,Y,Z] = meshgrid(x,x,x);

D = sqrt( X.^2 + Y.^2 ); % distance to central axis
a = 0.5; b = 0.2;
R = a + b*cos(pi*Z);

M = double( D<=R ); %  + 0.1*(1-2*rand(n,n,n));
% W = compute_edge_energy(M, 3);
W = abs(D-R)+0.1;
W = 1./W;

plot_volumetric_data(W);

% start/end points
zstart = 0.9;
rstart = a + b*cos(pi*zstart);
zend = -0.9;
rend = a + b*cos(pi*zend);

nstart = 100;
theta = 0:2*pi/nstart:2*pi-2*pi/nstart;
start_points = [rend*cos(theta); rend*sin(theta); zstart*ones(1,nstart)];
start_points = 0.5+0.5*start_points;
start_points = round(1+start_points*(n-1));

nend = 12*4;
theta = 0:2*pi/nend:2*pi-2*pi/nend;
end_points = [rend*cos(theta); rend*sin(theta); zend*ones(1,nend)];
end_points = 0.5+0.5*end_points;
end_points = round(1+end_points*(n-1));

% find unique points
id = start_points(1,:)+start_points(2,:)*pi+start_points(3,:)*pi^2;
[tmp,I] = unique(id);
start_points = start_points(:,I);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using line to point distance function
% perform marching
disp('Performing front propagation.');
clear options;
options.nb_iter_max = Inf;
tic;
[D,S] = perform_fast_marching(W, start_points, options);
toc;

disp('Extracting path.');
paths1 = compute_geodesic(D,end_points, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using point to point distance function
% perform marching
start_points = [end_points(1,:); end_points(2,:); ones(1,nend)*start_points(3,1)];
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
plot_fast_marching_3d(W,[],paths1,[],[]);
colormap gray(256);

subplot(1,2,2);
plot_fast_marching_3d(W,[],paths2,[],[]);
colormap gray(256);


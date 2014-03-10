% test for constrainted path planing.
map_name = 'data/room1.png';        % a small map
map_name = 'data/room.png';         % warning : this is a *big* map, expect large computation time.

M = imread(map_name);
M = 1-rescale( double(M) );
n = size(M,1);

% number of direction sampled
ntheta = 21;
% size of the object
obj = [0.2,0.03];

M1 = generate_constrained_map(M,ntheta,obj,pi);
W = 1-M1+0.001;

% pick start/end/center points
[start_points,end_points] = pick_start_end_point(M);
start_theta = floor(ntheta*rand)+1;      % (ntheta+1)/2
end_theta = floor(ntheta*rand)+1;        % (ntheta+1)/2
start_points = [start_points;start_theta];
end_points = [end_points;end_theta];


clear options;
options.nb_iter_max = Inf;
options.end_points = end_points;
disp('Performing 3D front propagation.');
[D,S] = perform_fast_marching(W, start_points, options);
options.Tmax = 3*n;
disp('Extracting shortest path.');
path = compute_geodesic(D,end_points, options);

path(:,3) = (path(:,3)-1)/ntheta*pi;

clf;
plot_constrained_path_planing(M,path,obj,40);
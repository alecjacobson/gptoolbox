% test for constrainted path planing using an heuristic.
%

save_image = 0;

name = 'room_big';         % warning : this is a *big* map, expect large computation time.
name = 'room';         % warning : this is a *big* map, expect large computation time.
name = 'room1';        % a small map
map_name = ['data/' name '.png'];

rep = 'images/';
rep_eps = 'images/eps/';

M = imread(map_name);
M = 1-rescale( double(M) );
n = size(M,1);

% number of direction sampled
ntheta = 41;
% size of the object
obj = [0.2,0.03];

M1 = generate_constrained_map(M,ntheta,obj,pi);
W = 1-M1+0.001;

% pick start/end/center points
if ~exist('start_points')
    [start_points,end_points] = pick_start_end_point(M);
    start_theta = floor(ntheta*rand)+1;      % (ntheta+1)/2
    end_theta = floor(ntheta*rand)+1;        % (ntheta+1)/2
    start_theta = 1;
    end_theta = 1;
    start_points = [start_points;start_theta];
    end_points = [end_points;end_theta];
end

clear options;
options.nb_iter_max = Inf;
options.end_points = end_points;
disp('Performing 3D front propagation.');
[D,S] = perform_fast_marching(W, start_points, options);
disp('Extracting shortest path.');
path = compute_geodesic(D,end_points);
path1(:,3) = (path(:,3)-1)/ntheta*pi;

clf;
plot_constrained_path_planing(M,path1,obj,40);


test_type = 'weight';
if strcmp(test_type,'weight')
    weight_list = [0 0.3 0.5 0.7 0.8 0.9 1 1.1];
    reduc_factor_list = 1 + weight_list*0;
else
    reduc_factor_list = 0:0.05:1.1;
    weight_list = 0.5 + reduc_factor_list*0;
end




clf;
nbr_paths = length(weight_list);
for i = 1:nbr_paths
    options.reduc_factor = reduc_factor_list(i);
    options.weight = weight_list(i);
    options.start_points = start_points;
    disp( sprintf('Performing FM*, weight=%d%%, resol=%d%%.', round(100*options.weight), round(100*options.reduc_factor)) );
    [D,S] = perform_fmstar_3d(W, start_points,end_points, options);    
    % extract path
    disp('Extracting Paths');
    path = compute_geodesic(D,end_points);
    clf;
    options.plot_planes = 0;
    plot_fast_marching_3d(W,S,path,start_points,end_points, options );
    if save_image
        str = [rep name '_front_' num2str( round(weight_list(i)*100) )];
        saveas( gcf, [str '.png'], 'png' );
    end
end
% test for the heuristically driven 
%   Fast Marching algorithm, 
%   aka FM* with landmark points
%
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% speed fonction W

rep = 'images/landmark_fmstar/';
if exist(rep)~=7
    mkdir(rep);
end
rep_eps = [rep 'eps/'];
if exist(rep_eps)~=7
    mkdir(rep_eps);
end

n = 300;

save_image = 1;
test_type = 'landmarks';
clear options;

name = 'stephanodiscusniagarae';
name = 'map';
name = 'constant';
name = 'bump';
name = 'gaussian';
name = 'cavern';
name = 'road2';
name = 'mountain';

[M,W] = load_potential_map(name, n);


% pick start/end/center points
nbr_start_points = 25;
if ~exist('start_points') || length(start_points)~=nbr_start_points
    for i=1:nbr_start_points
        [start_points{i},end_points{i}] = pick_start_end_point(M);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform FM* algorithm
landmark_init = 'farthest';
landmark_init = 'rand';

nbr_landmarks_list = [1:1:10, 12:2:20, 25:5:40];
weight = 1;

paths = {};
paths_old = {};
fronts = {};

nbr_paths = length(nbr_landmarks_list);

options.point_size = 20;
options.end_point_style = 'b.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extraction without heuristic
front_size = zeros(nbr_paths+1, 1);
for s=1:nbr_start_points
    options.end_points = end_points{s};
    [D,S] = perform_fast_marching(W, start_points{s}, options);
    path_0 = compute_geodesic(D,end_points{s},options);
    plot_fast_marching_2d(M,S,path_0, start_points{s}, end_points{s}, options);
    if s<=9
        str = [name '_' test_type  '_00lanmarks_test' num2str(s)];
        saveas(gcf, [rep_eps str '.eps'], 'eps');
        saveas(gcf, [rep str '.png'], 'png');
    end
    front_size(1) = front_size(1) + sum( S(:)<=0 );
end
front_size(1) = front_size(1)/nbr_start_points;


clf;
S_accum = zeros(n);
landmark = []; DL = [];
for i = 1:nbr_paths
    nbr_landmarks = nbr_landmarks_list(i);
    % update the lanmarks point list
    prev_nbr_landmarks = size(landmark,2);
    m = nbr_landmarks - prev_nbr_landmarks;
    landmark = [landmark, floor( rand(2, m)*(n-1) )+1];
    % update the distance maps
    DL2 = DL; DL = zeros(n,n,nbr_landmarks);
    if prev_nbr_landmarks>=1
        DL(:,:,1:prev_nbr_landmarks) = DL2;
    end
    for k = prev_nbr_landmarks+1:nbr_landmarks
        DL(:,:,k) = perform_fast_marching(W, landmark(:,k));
    end
    % perform the FM*
    for s=1:nbr_start_points
        options.weight = weight;
        options.start_points = start_points{s};
        options.heuristic = compute_heuristic_landmark(DL,end_points{s});
        disp( sprintf('Performing FM*, weight=%d%%, nb.landmarks=%d.', round(100*weight), nbr_landmarks ) );
        [D,S] = perform_fmstar_2d(W, start_points{s},end_points{s}, options);
        % size of the front
        front_size(i+1) = front_size(i+1) + sum( S(:)<=0 );
        if s==1
            % extract path
            path = compute_geodesic(D,end_points{s},options);
            num_landmark_str = num2str(nbr_landmarks);
            if nbr_landmarks<10
                num_landmark_str = ['0' num_landmark_str];
            end
            if 1 % i<=5
                plot_fast_marching_2d(M,S,path, start_points{s}, [end_points{s},landmark], options);
                str = [name '_' test_type  '_'  num_landmark_str 'landmarks'];
                saveas(gcf, [rep_eps str '.eps'], 'eps');
                saveas(gcf, [rep str '.png'], 'png');
            end
        end
    end
    front_size(i+1) = front_size(i+1)/nbr_start_points;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extraction with exact heuristic
options.heuristic = perform_fast_marching(W, end_points{1});
options.start_points = start_points{1};
options.end_points = end_points{1};
[D,S] = perform_fmstar_2d(W, start_points{1},end_points{1}, options);
path_inf = compute_geodesic(D,end_points{1},options);
plot_fast_marching_2d(M,S,path_inf, start_points{1}, end_points{1}, options);
str = [name '_' test_type  '_INFlanmarks'];
saveas(gcf, [rep_eps str '.eps'], 'eps');
saveas(gcf, [rep str '.png'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the speed up saving
% plot error
x = [0, nbr_landmarks_list];
y = 100-front_size/front_size(1)*100;
sel = find(x>1 & x<=25);
plot(x(sel), y(sel));
axis( [min(x(sel)) max(x(sel)), 70, 98] );
xlabel('Nbr.landmarks');
ylabel('Computation saving (in %)');
if save_image
    str = [name test_type '_computation_savings'];
    saveas(gcf, [rep_eps str '.eps'], 'eps');
    saveas(gcf, [rep str '.png'], 'png');
end
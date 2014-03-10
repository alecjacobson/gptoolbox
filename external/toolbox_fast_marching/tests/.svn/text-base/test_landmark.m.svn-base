% test for landmark heuristic computation

n = 300;

rep = 'images/landmark_distance/';
if exist(rep)~=7
    mkdir(rep);
end
rep_eps = [rep 'eps/'];
if exist(rep_eps)~=7
    mkdir(rep_eps);
end

name = 'bump';
name = 'map';
name = 'mountain';
name = 'stephanodiscusniagarae';
name = 'cavern';
name = 'gaussian';
name = 'road2';
name = 'constant';

save_images = 0;

nbr_landmarks = 10;
landmark_init = 'rand';
landmark_init = 'farthest';

[M,W] = load_potential_map(name, n);

% just a string representing the number of landmarks points
num_str = num2string_fixeddigit(nbr_landmarks, 3);
marker_size = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick starting point
if ~exist('start_points')
    [start_points,end_points] = pick_start_end_point(M);
end

disp('Performing front propagation.');
[D,S] = perform_fast_marching(W, start_points);

% compute landmark
disp('Initializing Landmarks');
switch lower(landmark_init)
    case 'rand'
        landmark = floor( rand(2, nbr_landmarks)*(n-1) )+1;
    case 'farthest'
        landmark = perform_farthest_point_sampling( W, [], nbr_landmarks );
    case 'farthestboundary'
        landmark = perform_farthest_point_sampling_boundary( W, landmark, 1 );
end

% compute distance to landmark
DL = zeros(n,n,nbr_landmarks);
for i=1:nbr_landmarks
    disp('Performing front propagation.');
    [DL(:,:,i),S] = perform_fast_marching(W, landmark(:,i));
end

% compute the approximated distance
[D1,Z] = compute_heuristic_landmark(DL,start_points);

clf;
if ~save_images
    subplot(2,2,1);
end
hold on;
imagesc(W');
plot_scattered(landmark);
axis tight; axis square; axis off;
hold off;
title('Potential with landmarks');
colormap jet(256);
if save_images
    axis image; axis off;
    str = [name '_' num_str 'landmarks_' landmark_init '_potential'];
    saveas(gcf, [rep str '.png'], 'png');
    saveas(gcf, [rep_eps str '.eps'], 'eps');
end

if ~save_images
    subplot(2,2,3);
end
hold on;
imagesc(D');
contour(D', 'k');
axis tight; axis square; axis off;
plot(start_points(1,:), start_points(2,:), 'kx', 'MarkerSize', 10);
hold off;
title('Distance to point');
colormap jet(256);
if save_images
    axis image; axis off;
    str = [name '_' num_str 'landmarks_' landmark_init '_realdist'];
    saveas(gcf, [rep str '.png'], 'png');
    saveas(gcf, [rep_eps str '.eps'], 'eps');
end

if ~save_images
    subplot(2,2,4);
end
hold on;
imagesc(D1');
contour(D1', 'k');
axis tight; axis square; axis off;
plot(start_points(1,:), start_points(2,:), 'kx', 'MarkerSize', 10);
plot(landmark(1,:), landmark(2,:), 'k.', 'MarkerSize', marker_size);
hold off;
title('Approximate distance to point');
colormap jet(256);
if save_images
    axis image; axis off;
    str = [name '_' num_str 'landmarks_' landmark_init '_approxdist'];
    saveas(gcf, [rep str '.png'], 'png');
    saveas(gcf, [rep_eps str '.eps'], 'eps');
end

% str = [rep name '_' num2str(nbr_landmarks) 'landmarks_' landmark_init];
% saveas(gcf, [str '.png'], 'png');

return;

% plot the influence zones
clf;
hold on;
imagesc(Z');
plot_scattered(start_points);
axis off;
hold off;
title('Approximate distance to point');
colormap jet(256);


    
return;


%%% ANIMATION! %%%
r = 0.3;
nb = 200;
i = 0;
for theta = (0:1/nb:1-1/nb)*2*pi
    i = i+1;
    start_points = 0.5 + r * [cos(theta);sin(theta)];
    start_points = round( start_points*(n-1) )+1;
    [D1,Z] = compute_distance_landmark(DL,start_points);
    clf;
    imagesc(Z);
    axis image;
    axis off;
    colormap jet(256);
    
    num_str = num2str(i);
    if i<10
        num_str = ['0' num_str];
    end
    if i<100
        num_str = ['0' num_str];
    end
    str = ['anim/' name '_' num_str];
    saveas(gcf, [str '.png'], 'png');
end
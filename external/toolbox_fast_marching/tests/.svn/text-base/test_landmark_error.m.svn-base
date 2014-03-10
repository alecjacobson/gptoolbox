% test for influence of sampling strategy on the 
% decreasing of the error


n = 100;

rep = 'images/landmark_error/';
rep_eps = [rep  'eps/'];
if exist(rep)~=7
    mkdir(rep);
end
if exist(rep_eps)~=7
    mkdir(rep_eps);
end

marker_size = 20;
contour_width = 2;
nb_contours = 12;
save_image = 1;
save_image_eps = 0;

if ~exist('name')
    name = 'bump';
    name = 'map';
    name = 'stephanodiscusniagarae';
    name = 'cavern';
    name = 'gaussian';
    name = 'constant';
    name = 'road2';
    name = 'mountain';
end

[M,W] = load_potential_map(name, n);

switch lower(name)
    case 'road2'
        A = [0.0352    0.7638; 0.5377    0.4271];
    case 'mountain'
        A = [0.1005    0.8442; 0.6332    0.3116];
    case 'constant'
        A = [0.1 0.1; 0.7 0.9];
end

% get two samples points
if exist('A')
    start_points = round(A(:,1)*(n-1)+1);
    end_points = round(A(:,2)*(n-1)+1);
else
    [start_points,end_points] = pick_start_end_point(M);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save real distance function
options.end_points = end_points;
[D,S] = perform_fast_marching(W, start_points, options);
path = compute_geodesic(D,end_points);
clf;
plot_fast_marching_2d(M,S,path, start_points, end_points);
% save to file
if save_image
    str = [name '_original_path'];
    saveas(gcf, [rep str '.png'], 'png');
    if save_image_eps
        saveas(gcf, [rep_eps str '.eps'], 'epsc');
    end
end

% plot distance function
[D,S] = perform_fast_marching(W, end_points);
clf;
hold on;
m = sort(D(:)); m = m(end-10);
imagesc(min(D,m)');
contour(min(D,m)', nb_contours, 'k', 'LineWidth', contour_width);
axis tight; axis image; axis off;
plot(start_points(1,:), start_points(2,:), 'kx', 'MarkerSize', 10);
hold off;
colormap jet(256);
axis image; axis off;
% save to file
if save_image
    str = [name '_original_distance'];
    saveas(gcf, [rep str '.png'], 'png');
    if save_image_eps
        saveas(gcf, [rep_eps str '.eps'], 'epsc');
    end
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% base point for error computation %%%%%%
nbr_base_points = 200;
% base point position
base_points = floor( rand(2, nbr_base_points)*n ) + 1;
% record distances to landmarks
DL_base = zeros(n,n,nbr_base_points);  
% compute real distance to these base ponts
fprintf('Compute distance to base points ');
for s=1:nbr_base_points
    fprintf('.');
    [DL_base(:,:,s),Z] = perform_fast_marching(W, base_points(:,s));
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% error decreasing %%%%%%
% landmark_init_list = {'farthestboundary','rand', 'farthest','trialerror'}; 
landmark_method_list = {'proxy','align'};
landmark_init_list{1} = {'rand'}; 
landmark_init_list{1} = {}; 
landmark_init_list{2} = {'rand', 'farthestboundary'}; % 'trialerror' 'farthestboundary'
err_list = {};
% for recording error
nbr_landmark_list{1} = [5:10:128, 128:32:640];
nbr_landmark_list{2} = [1, 2:2:32];
% for recording images
display_list{1} = [65, 128, 256, 384, 512, 640];
display_list{2} = [1, 2, 4, 8, 16, 32];

% performing test
err_list = {};
err_list_abscice = {};
lgd = {};
k = 0;
for k1=1:length(landmark_method_list)
for k2=1:length(landmark_init_list{k1})
    landmark_method = landmark_method_list{k1};
    landmark_init = landmark_init_list{k1}{k2};
    lgd{end+1} = [landmark_method '+' landmark_init];
    
    err = [];
    % landmark poisition
    landmark = [];
    % record distances to landmarks
    DL_landmark = [];
    
    for nbr_landmarks=nbr_landmark_list{k1}
        
        disp( ['Testing method with ' num2str(nbr_landmarks) ' landmarks'] );
        while size(landmark,2)<nbr_landmarks
            % update landmarks points
            switch lower(landmark_init)
                case 'rand'
                    landmark = [landmark, floor( rand(2, 1)*(n-1) )+1];
                case 'farthest'
                    landmark = perform_farthest_point_sampling( W, landmark, 1 );
                case 'trialerror'
                    landmark = perform_farthest_landmark_sampling( W, landmark, DL_landmark, base_points, DL_base, 1 );
                case 'farthestboundary'
                    landmark = perform_farthest_point_sampling_boundary( W, landmark, 1 );
                case 'farthestunif'
                    landmark = perform_farthest_point_sampling_boundary( W, landmark, 1, 'unif' );
            end
            
            % update distance map
            [DL1,S] = perform_fast_marching(W, landmark(:,end));
            DL_landmark = cat(3,DL_landmark,DL1);
        end
        
        % compute a mean error
        e = 0;
        for s=1:nbr_base_points
            % distance using landmarks
            D1 = compute_distance_landmark(base_points(:,s), DL_landmark, landmark, landmark_method);
            % real distance
            D = DL_base(:,:,s);
            e = e + mean( (D1(:)-D(:)).^2 );
        end        
        err = [err; e / nbr_base_points];
        
        % save images
        if ~isempty(find(display_list{k1}==nbr_landmarks))
            
            % base name for file saving 
            axis image; axis off;
                
            % save distance function
            D1 = compute_distance_landmark(end_points, DL_landmark, landmark, landmark_method);
            clf;
            hold on;
            m = sort(D1(:)); m = m(end-10);
            imagesc(min(D1,m)');
            if strcmp(landmark_method, 'align')
                contour(min(D1,m)', nb_contours, 'k','LineWidth', contour_width);
            end
            axis tight; axis image; axis off;
            plot(end_points(1,:), end_points(2,:), 'kx', 'MarkerSize', 10);
            % plot(landmark(1,:), landmark(2,:), 'k.', 'MarkerSize', marker_size);
            hold off;
            colormap jet(256);
            % save to file
            if save_image
                num_str = num2string_fixeddigit(nbr_landmarks, 3);
                str = [name '_' landmark_init '_' landmark_method '_distance_' num_str 'ldm'];
                saveas(gcf, [rep str '.png'], 'png');
                if save_image_eps
                    saveas(gcf, [rep_eps str '.eps'], 'epsc');
                end
            end
            
            % save geodesic extraction + saving
            options.heuristic = D1;
            options.weight = 1;
            [D,S] = perform_fmstar_2d(W, start_points,end_points, options);
            path = compute_geodesic(D,end_points);
            clf;
            options.point_size = {12, marker_size};
            plot_fast_marching_2d(M,S,path,start_points,[end_points,landmark], options);
            % save to file
            if save_image
                str = [name '_' landmark_init '_' landmark_method '_path_' num_str 'ldm'];
                saveas(gcf, [rep str '.png'], 'png');
                if save_image_eps
                    saveas(gcf, [rep_eps str '.eps'], 'epsc');
                end
            end
        end
    end
    err_list{end+1} = err;
    err_list_abscice{end+1} = nbr_landmark_list{k1};
end
end

% plot error
clf;
hold on;
for i=1:length(err_list)
    x = err_list_abscice{i};
    if max(x)==max(nbr_landmark_list{1})
        x = x.^2+n^2;
    else
        x = x*n^2;
    end
    plot(log2(x), log2(err_list{i}), get_color_from_index(i) );
end
hold off;
legend(lgd);
xlabel('log_{10}(#landmarks)');
ylabel('log_{10}(Error)');
axis tight;
if save_image
    str = [name '_error'];
    saveas(gcf, [rep str '.png'], 'png');
    if save_image_eps
        saveas(gcf, [rep_eps  str '.eps'], 'epsc');
    end
end
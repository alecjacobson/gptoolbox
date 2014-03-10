% test for influence of sampling strategy on the 
% decreasing of the error

n = 80;

rep = 'images/distance_approx/';
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
    name = 'road2';
    name = 'mountain';
    name = 'constant';
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
D = my_eval_distance(W, end_points);
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
nbr_base_points = 50;
% base point position
base_points = floor( rand(2, nbr_base_points)*n ) + 1;
% record distances to landmarks
DL_base = zeros(n,n,nbr_base_points);  
% compute real distance to these base ponts
fprintf('Compute distance to base points ');
for s=1:nbr_base_points
    fprintf('.');
    DL_base(:,:,s) = my_eval_distance(W, base_points(:,s));
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% error decreasing %%%%%%
err_list = {};
% for recording error
nbr_landmark_list = [16:16:128];
% for recording images
display_list = [16 32 64 128 256];
landmark_init_list = {'rand', 'farthestboundary'}; 
landmark_init_list = {'rand'}; % 'farthest' 'trialerror' 'farthestboundary'

% performing test
err_list = {};
err_list_abscice = {};
lgd = {};
k = 0;

ldm_order_list = 0:5;
err = zeros(length(ldm_order_list), length(nbr_landmark_list));

for k=1:length(landmark_init_list)
    
    landmark_init = landmark_init_list{k};
    
    % landmark position
    landmark = [];
    % record distances to landmarks
    DL_landmark = [];
    
    for t=1:length(nbr_landmark_list)
        nbr_landmarks = nbr_landmark_list(t);
        
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
            DL1 = my_eval_distance(W, landmark(:,end));
            DL_landmark = cat(3,DL_landmark,DL1);
        end
        
        % compute a mean error
        for s=1:length(ldm_order_list)
            ldm_order = ldm_order_list(s);
            landmark_method = ['multiproxy' num2str(ldm_order)];
            e = 0;
            for pt=1:nbr_base_points
                % distance using landmarks
                D1 = compute_distance_landmark(base_points(:,pt), DL_landmark, landmark, landmark_method);
                % real distance
                D = DL_base(:,:,pt);
                e = e + sqrt( mean( (D1(:)-D(:)).^2 ) );
            end
            err(s,t) = e/nbr_base_points;
        end
        
        % save images
        if ~isempty(find(display_list==nbr_landmarks))
            for s=1:length(ldm_order_list)
                ldm_order = ldm_order_list(s);
                landmark_method = ['multiproxy' num2str(ldm_order)];
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
                    str = [name '_' landmark_init '_'  num_str 'ldm_' landmark_method '_distance'];
                    saveas(gcf, [rep str '.png'], 'png');
                    if save_image_eps
                        saveas(gcf, [rep_eps str '.eps'], 'epsc');
                    end
                end
            end
        end
    end
end

% plot error
clf;
hold on;
for i=1:size(err,1)
    x = nbr_landmark_list;
%    x = x.^2+i*n^2;
    plot(log10(x), log10(err(i,:)), get_color_from_index(i) );
end
hold off;
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

% plot error
clf;
hold on;
for i=1:size(err,1)
    x = nbr_landmark_list;
    if ldm_order_list(i)==0
        x = x*n^2;  % align method
    else
        x = x.^2+ldm_order_list(i)*n^2;
    end
    plot(log10(x), log10(err(i,:)), get_color_from_index(i) );
end
hold off;
xlabel('log_{10}(memory)');
ylabel('log_{10}(Error)');
axis tight;
if save_image
    str = [name '_error_memory'];
    saveas(gcf, [rep str '.png'], 'png');
    if save_image_eps
        saveas(gcf, [rep_eps  str '.eps'], 'epsc');
    end
end
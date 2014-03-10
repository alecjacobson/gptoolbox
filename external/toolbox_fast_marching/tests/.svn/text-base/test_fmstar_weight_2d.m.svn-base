% test for the heuristically driven 
%   Fast Marching algorithm, 
%   aka FM*
%
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% speed fonction W

n = 500;


rep = 'images/multires_weight/';
if exist(rep)~=7
    mkdir(rep);
end
rep_eps = [rep 'eps/'];
if exist(rep_eps)~=7
    mkdir(rep_eps);
end


save_image = 1;
test_type = 'resolution';
test_type = 'weight';

name = 'bump';
name = 'stephanodiscusniagarae';
name = 'map';
name = 'cavern';
name = 'road2';
name = 'constant';
name = 'gaussian';
name = 'mountain';

[M,W] = load_potential_map(name, n);

% pick start/end/center points
if ~exist('start_points')
    [start_points,end_points] = pick_start_end_point(M);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform FM* algorithm
if strcmp(test_type,'weight')
    weight_list = 0:0.05:1; % 1.1;
    reduc_factor_list = 1 + weight_list*0;
else
    reduc_factor_list = 0.05:0.05:1;
    weight_list = 0.5 + reduc_factor_list*0;
end

paths = {};
paths_old = {};
fronts = {};

nbr_paths = length(weight_list);

clf;
S_accum = zeros(n);
front_size = zeros(nbr_paths, 1);
for i = 1:nbr_paths
    options.reduc_factor = reduc_factor_list(i);
    options.weight = weight_list(i);
    options.start_points = start_points;
    disp( sprintf('Performing FM*, weight=%d%%, resol=%d%%.', round(100*options.weight), round(100*options.reduc_factor)) );
    [D,S] = perform_fmstar_2d(W, start_points,end_points, options);    
    S_accum = S_accum + (S+1)/2;
    % size of the front
    front_size(i) = sum( S(:)<=0 );
    % extract path
    disp('Extracting Paths');
    paths{i} = compute_geodesic(D,end_points,options);
    if i==1
        plot_fast_marching_2d(M,S,paths{i}, start_points, end_points, options);
    end
    paths{i} = (paths{i}-1)/n;
    % extract front
    c_list = perform_curve_extraction(S,0);
    fronts{i} = c_list{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the level sets
clf;
hold on;
x = linspace(0,1,n);
imagesc(x,x,M');
axis image;
axis off;
plot_curve(fronts, [], 'r', 1);
plot_curve(paths{end}, [], 'b', 2);
hold off;
colormap gray(256);
if 0 % save_image
    str = [name '_' test_type  '_fronts'];
    saveas(gcf, [rep str '.png'], 'png');
    saveas(gcf, [rep_eps str '.eps'], 'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the geodesics
clf;
hold on;
x = linspace(0,1,n);
imagesc(x,x,M');
axis image;
axis off;
plot_curve({paths{2:end}}, [], 'k', 1);
plot_curve(paths{1}, [], 'b', 2);
hold off;
colormap gray(256);
if save_image
    str = [name '_' test_type '_paths'];
    saveas(gcf, [rep str '.png'], 'png');
    saveas(gcf, [rep_eps str '.eps'], 'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the hausdorf distance
err_rms = zeros(nbr_paths, 1);
err_inf = zeros(nbr_paths, 1);
for i=1:nbr_paths
    [err_rms(i), err_inf(i)] = compute_hausdorff_distance(paths{i}, paths{1});
end
% normalize the RMS in %
err = {};
err{1} = 100 * err_rms / norm( start_points - end_points, 'fro' );
err{2} = 100 * err_inf / norm( start_points - end_points, 'fro' );
str1 = {'rms', 'inf'};
% plot error
for i=1:2
    if strcmp(test_type,'weight')
        weight_max = 1;
        I = find(weight_list<=weight_max);
        plot( weight_list(I)*100, err{i}(I) );
        axis tight;
        xlabel('Heuristic (in %)');
        ylabel( ['Hausdorff ' upper(str1{i}) ' (in % B.B.diag.)'] );
    else
        plot( reduc_factor_list*100, err{i} );
        axis tight;
        xlabel('% reduction');
        ylabel( ['Hausdorff ' upper(str1{i}) ' (in % B.B.diag.)'] );
    end
    if save_image
        str = [name '_' test_type '_error_' str1{i}];
        saveas(gcf, [rep_eps str '.eps'], 'epsc');
        saveas(gcf, [rep str '.png'], 'png');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the memory saving
% plot error
if strcmp(test_type,'weight')
    plot(weight_list*100, 100-front_size/front_size(1)*100);
    axis( [0, max(weight_list*100), 0, 100] );
    xlabel('Heuristic (in %)');
    ylabel('Memory saving (in %)');
    if save_image
        str = [name '_' test_type '_memory'];
        saveas(gcf, [rep_eps str '.eps'], 'eps');
        saveas(gcf, [rep str '.png'], 'png');
    end
end
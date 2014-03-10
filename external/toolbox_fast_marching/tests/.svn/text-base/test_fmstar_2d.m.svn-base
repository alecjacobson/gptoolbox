% test for the heuristically driven 
%   Fast Marching algorithm, 
%   aka FM*
%
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% speed fonction W

save_image = 1;

n = 600;

rep = 'images/multires/';
if exist(rep)~=7
    mkdir(rep);
end
rep_eps = [rep 'eps/'];
if exist(rep_eps)~=7
    mkdir(rep_eps);
end

% type = 'constant';
energy_type = 'image';
energy_type = 'edge';

propagation_type = 'normal';
propagation_type = 'circular';

geodesic_type = 'group';
geodesic_type = 'punctual';
group_witdh = 40;    % in pixels

name = 'bump';
name = 'cavern';
name = 'map';
name = 'mountain';
name = 'gaussian';
name = 'cavern';
name = 'road2';
name = 'stephanodiscusniagarae';
name = 'diatom2';

[M,W] = load_potential_map(name, n);
    

% pick start/end/center points
if strcmp(propagation_type, 'normal')
    [start_point,end_point] = pick_start_end_point(M);
else
    [start_point,end_point,center_point] = pick_start_end_point(M);
    options.center_point = center_point;
end

if strcmp(geodesic_type, 'group')
    start_points = repmat(start_point, 1, 2*group_witdh+1);
    start_points = start_points + [-group_witdh:group_witdh; zeros(1,2*group_witdh+1)];

    end_points = repmat(end_point, 1, 2*group_witdh+1);
    end_points = end_points + [-group_witdh:group_witdh; zeros(1,2*group_witdh+1)]; 
else
    end_points = end_point;
    start_points = start_point;
end

s = 3;
epsi = 0.000001;
if strcmp(energy_type, 'edge')
    if strcmp(propagation_type, 'normal')
        W = compute_edge_energy(M, 5);
    else
        W = compute_edge_energy(M, s, epsi, center_point);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform FM* algorithm
weight_list = [0,0.5,0.7,0.8,0.9];
weight_list = [0,0.5,0.8,0.9,0.95, 1, 1.02, 1.05];
reduc_factor = 0.8;
k = 0;

options.propagation_type = propagation_type;
options.Tmax = sum(size(M))*1.2;
options.start_points = start_points;

clf;
for weight = weight_list
    k = k+1;
    disp('Performing FM*');
    options.reduc_factor = reduc_factor;
    options.weight = weight;
    [D,S] = perform_fmstar_2d(W, start_points,end_points, options);    
    disp('Extracting Paths');
    
    % find the end point with minimum action
    d = [];
    for i=1:size(end_points,2)
        d = [d D(end_points(1,i),end_points(2,i))];
    end
    [tmp,I] = min(d);
    end_point_min = end_points(:,I(1));
    
    opt.stepsize = 0.2;
    path = compute_geodesic(D,end_point_min, opt);
    % subplot(1, length(weight_list), k);
    clf;
    plot_fast_marching_2d(M,S,path,start_point,end_point_min);
    % colormap gray(256);
    
    if save_image
        numstr = sprintf('%.2f', weight);
        str = ['fmstar_2d_', name, '_', propagation_type, '_', geodesic_type, '_', numstr];
        saveas(gcf, [rep_eps str '.eps'], 'epsc');
        saveas(gcf, [rep str '.png'], 'png');
        if 0
        warning off;
        imwrite(rescale(S), [rep 'STATE_' str '.png'], 'png' );
        warning on;
        end
    end
end


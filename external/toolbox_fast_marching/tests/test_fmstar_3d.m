% test for the heuristically driven 
%   Fast Marching algorithm, 
%   aka FM* in 3D
%
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% speed fonction W

save_image = 0;

n = 100;

type = 'heart';
type = 'constant';
type = 'gaussian';

reduc_factor = 0.5;

if strcmp(type, 'constant')    
    W = ones(n,n,n);
    W(1) = 0.99;
    M = W;
elseif strcmp(type, 'gaussian')
    x = -1:2/(n-1):1;
    [X,Y,Z] = meshgrid(x,x,x);
    sigma = 0.8;
    W = 1./(1 + exp( -(X.^2+Y.^2+Z.^2)/sigma^2 ) );
    M = W;
elseif strcmp(type, 'heart')
    load('heart1');
    a = [60, 40, 20];
    b = [180,160, 120];
    M = heart1(a(1):b(1),a(2):b(2),a(3):b(3));
    clear heart1;
    M = perform_image_resize(M, n,n,n );
    epsi = 0.05; s = 3;
    W = compute_edge_energy(M,s,epsi);
end


n = length(W);

k = 5;
AZ = 125;
weight_list = [0,0.5,1,1.1];
if strcmp(type, 'constant')
    start_points = [n/2;n/2;n/2];
    end_points = [k;n-k;n-k];
    weight_list = [0,0.5,1,1.2];
elseif strcmp(type, 'gaussian')
    start_points = [n-k; k; k];
    end_points = [k;n-k;n-k];
    AZ = 125+180;
    weight_list = [0,0.5,1,1.1,1.2];
elseif strcmp(type, 'heart')
    AZ = 156;
    start_points = [52; 66; 72];
    end_points = [53; 53; 8];
end

% plot the data
plot_volumetric_data( M );
view(AZ, 30);
colormap jet(256);
if save_image
    saveas(gcf, ['fmstar_3d_',type, '_map'], 'png');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform FM* algorithm

reduc_factor = 1;
clf;
k = 0;
for weight = weight_list
    k = k+1;
    disp('Performing FM*');
    [D,S] = perform_fmstar_3d(W, start_points,end_points, reduc_factor,weight);    
   
    path = compute_geodesic(D,end_points);
    
    options.plot_planes = 0;
    clf;
    plot_fast_marching_3d(M,S,path, start_points, end_points, options);
    colormap gray(256);
    axis([1 n 1 n 1 n]);    
    camlight;
    view(AZ, 30);
    if strcmp(type, 'heart')
        colormap jet(256);
    end
    
    if save_image
        saveas(gcf, ['fmstar_3d_',type, '_', num2str(k)], 'png');
    end
end


% saveas(gcf, ['fmstar_2d_',type], 'eps');
% saveas(gcf, ['fmstar_2d_',type], 'png');
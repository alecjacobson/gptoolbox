% test for anisotropic propagation in 2D

name = 'fixed-2d';
name = 'fixed-3d';
name = 'varying-2d';

rep = 'results/anisotropic-remeshing/';
if not(exist(rep))
    mkdir(rep);
end

% spacially varying 2D field
n = 200;
% create 2D vector field
s = [n n 1];
U = randn(n,n,2);
sigma = 30;
for it=1:10
    U = perform_vf_normalization( perform_blurring(U, sigma) );
end
U = perform_vf_normalization( U );

% test for various degree of anisotropy
% aniso_list = [0.01 0.05 .1 .2 .5 1];
aniso = .05;

% 3D field
V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
T = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );
D = perform_fast_marching(T, start_points);


clf;
hold on;
options.sub = round(n/15);
options.color = 'k';
plot_tensor_field(T1, D1, options);
h = plot(start_points(2,:),start_points(1,:), 'r.');
set(h, 'MarkerSize', 20);
colormap jet(256);

%    saveas(gcf, [rep name '-aniso-' num2str(ianiso) '.png'], 'png');

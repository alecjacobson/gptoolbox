% test for Fethalah code and Prados code
n = 500;

rep = 'results/anisotropic-feth/';
if not(exist(rep))
    mkdir(rep);
end

%% compute random tensor field
randn('seed', 12345);
U = randn(n,n,2);
sigma = (n/200)*30;
options.bound = 'per';
for it=1:10
    U = perform_vf_normalization( perform_blurring(U, sigma,options) );
end
U = perform_vf_normalization( U );

aniso_list = [1 .5 .2 .1 .05 .01 .001];

p = 5;
x = round( n/(p*2):n/p:n-n/(p*2) );
[Y,X] = meshgrid(x,x);
% start_points = round( start_points*(n-1)+1 );
start_points = cat(1, X(:)', Y(:)');    

s = randperm(size(start_points,2));
for i=1:length(aniso_list)

    aniso = aniso_list(i);
    V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
    T = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );
    Ti = perform_tensor_recomp(U,V, ones(n),ones(n)*1/aniso );

    options.use_feth_code = 0;
    tic;
    [D,S,Q] = perform_fast_marching(T, start_points, options);
    disp(['Prados: ' num2str(toc)]);

    hx = 1/n; hy = 1/n;
    tic;
    [D1,dUx,dUy, Vor, L] = fm2dAniso([hx;hy], Ti, start_points);
    disp(['Feth: ' num2str(toc)]);

    imageplot(s(Q), 'Prados', 1,2,1);
    imageplot(s(Vor+1), 'Feth', 1,2,2);
    colormap jet(256);
    saveas(gcf, [rep 'voronoi-' num2str(aniso) '.png'], 'png');
    
    clf;
    plot_tensor_field(T, perform_histogram_equalization(D, 'linear'));
    colormap jet(256);
    saveas(gcf, [rep 'distance-aniso-' num2str(aniso) '-prados.png'], 'png');
    
    clf;
    plot_tensor_field(T, perform_histogram_equalization(D1, 'linear'));
    colormap jet(256);
    saveas(gcf, [rep 'distance-' num2str(aniso) '-feth.png'], 'png');

end
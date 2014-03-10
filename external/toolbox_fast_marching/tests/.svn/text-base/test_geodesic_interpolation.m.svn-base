% test for data interpolation using geodesic distances

rep  = 'data/';
rep  = '';

name = 'mm';
name = 'cavern';
name = 'disk';
name = 'toto';
name = 'boat';
name = 'cartoon';

clear options;
options.null = 0;

%% compute a binary shape
n = 300;
if strcmp(name, 'toto')
    M0 = load_image([rep name],n);
    M = sum(M0,3);
    M = 1-(M==M(1));
    % compute the constraint map
    CM = zeros(n) + Inf; CM(M==0) = -Inf;
    options.constraint_map = CM; % constraint the propagation inside the shape
    W = ones(n);
    mode = 'interpolation';
else
    n0 = [];
    if strcmp(name, 'disk')
        n0 = n;
    end
    if strcmp(name, 'cartoon')
        n0 = 400;
    end
    M = load_image([rep name],n0);
    M = rescale(crop(M,n)); M0 = M;
    M = rescale(sum(M,3));
    W = compute_edge_energy(M, 1.2);
    W = 1./rescale(W,0.001,1);
    mode = 'colorization';
end

% RGB2COL matrix
CM = rand(3); CM(:,1) = 1;
[CM,R] = qr(CM); CM(:,1) = 1/3; CM = CM'; 

%% pick some points
b = 1;
points = [];
while b==1
    clf;
    hold on;
    imagesc(rescale(M)); axis image; axis off;
    if not(isempty(points))
        plot_scattered(points(end:-1:1,:));
    end
    colormap gray(256); axis ij;
    hold off;
    [y,x,b] = ginput(1);
    if b==1
        points(:,end+1) = [x;y];
    end
end
npoints = size(points,2);
% random colors (sqrt to make it more colorfull
f = rand(3,npoints);
f = ( f ./ repmat(sum(f,1),[3, 1]) ).^(1);
% use fixed palette
fc = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 0.2 0.7 1; 1 0.2 0.7; 0.7 1 0.2; 0.2 1 0.7]';
sel =  floor(rand(npoints,1)*6)+1;
f = fc(:,sel);

if strcmp(mode, 'colorization')
    % random CbCr
    f = CM*f; f = f(2:end,:);
end

%% perform interpolation
if strcmp(mode, 'colorization')
    sigma_list = 0.001/n;
    alpha_list = 3;
else
    % the smaller, the more interpolation is performed
    % the higher, the more diffusion it is
    sigma_list = [25 25 25 0.01 0.01 0.01]/n;
    % the higher, the more "voronoi"-like it is
    alpha_list = [1 3 50 1 3 50];
end

nsigma = length(sigma_list);
nrows = 1;
if nsigma>3
    nrows = 2;
end

clf;
for i=1:nsigma
    options.sigma = sigma_list(i);
    options.alpha = alpha_list(i);
    [A,G] = perform_geodesic_interpolation(W,points,f,options);
    
    if strcmp(mode, 'colorization')
        % cat the luminosity component
        A = cat(3, M,A);
        A = reshape(A, [n^2 3])';
        A = ( (CM^(-1)) * A )';
        A = reshape( A, [n n 3] );
    end
    A = clamp(A);
    subplot(nrows,ceil(nsigma/nrows), i);
    hold on;
    imagesc(rescale(A)); axis image; axis off; axis ij;
    plot_scattered(points(end:-1:1,:));
    hold off;
    title(['\alpha=' num2str(options.alpha) ', \sigma=' num2str(options.sigma)]);
end
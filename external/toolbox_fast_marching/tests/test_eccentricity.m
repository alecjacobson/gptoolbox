% compute the eccentricity of a 2D shape

path(path, 'toolbox/');
path(path, 'data/');

rep = 'images/';
name = 'bird';
name = 'apple';
name = 'device';
name = 'camel';
name = 'cavern';
name = 'chicken';
name = 'giraffe';

method = 'ecc';
method = 'mean';
metric = 'euclidean';
metric = 'geodesic';
options.metric = metric;
options.method = method;

% exponent for the computation
s = Inf;
if strcmp(method, 'mean')
    s = 1;
end

rep = 'results/eccentricity/';
if ~exist(rep)
    mkdir(rep);
end
testname = [method '-' metric];

n = 256;
M = rescale( load_image(name,n), 0,1 );
M = double(M>0.5);

% make sure pixels on the boundary are black
if M(1)==1
    M = 1-M;
end

% compare the geodesic distance to the euclidean distance
clf;
imagesc(M); axis image; axis off;
title('click on a point inside the shape');
% [y,x] = ginput(1);
x = 50; y = 50;
start_points = round([x y]');
W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;
options.constraint_map = [];
[D0,S,Q] = perform_fast_marching(W, start_points, options);
options.constraint_map = L;
[D1,S,Q] = perform_fast_marching(W, start_points, options);
D0(M==0) = 0; D1(M==0) = 0;

clf;
subplot(1,2,1);
hold on;
display_eccentricity(D0, 1);
h = plot(y,x, 'r.');
set(h, 'MarkerSize', 20);
axis ij;
title('Euclidean distance');
subplot(1,2,2);
hold on;
display_eccentricity(D1, 1);
h = plot(y,x, 'r.');
set(h, 'MarkerSize', 20);
axis ij;
title('Geodesic distance');

% saveas(gcf, [rep name '-geodesic-dist.eps'], 'epsc');


options.s = s;
options.nb_samples = 300; % subamples
[E,Ind] = compute_eccentricity_transform(M, options);
    
clf; 
subplot(1,2,1);
display_eccentricity(E, 1);
title('Eccentricity Transform');

if strcmp(method, 'ecc')
    [ep,I,J] = unique( Ind(:) ); ep(ep==Inf) = [];
    [eccx,eccy] = ind2sub([n n], ep);
    J = reshape(J,n,n);

    subplot(1,2,2);
    hold on;
    imagesc(J);
    h = plot(eccy,eccx,'r.');
    set(h, 'MarkerSize', 20);
    hold off;
    title('Eccentric points');
    axis image; axis off; axis ij;

    % saveas(gcf, [rep name '-' testname '.png'], 'png');
    % saveas(gcf, [rep name '-' testname '.eps'], 'epsc');

    % subplot(2,2,3);
    % display_eccentricity(E, 2);
end


% display the histograms
clf;
I = find(E>0);
v = rescale(E(I));
hist(v, 60);
title('Histogram');
axis tight;
% saveas(gcf, [rep name '-histogram.eps'], 'epsc');
saveas(gcf, [rep name '-histogram-' method '.png'], 'png');

%% save as pretty images
[Y,X] = meshgrid(1:n,1:n);
r = 0.02*n;
% save shape
warning off;
imwrite( 1-M, [rep name '-shape.png'], 'png' );
% save eccentricity
E0 = E; E0(M==0) = Inf;
E0 = convert_distance_color(E0);
% place a dot at the center
[tmp,i] = min(E(:) + (E(:)==0)*1e9); [x,y] = ind2sub([n n], i);
I = find( (X-x).^2 + (Y-y).^2 <= r^2 );
for i=1:3
    A = E0(:,:,i); A(I) = i==1;
    E0(:,:,i) = A;
    A = E0(:,:,i); A(I) = i==1;
    E0(:,:,i) = A;
end
imwrite( E0, [rep name '-shape-' method '.png'], 'png' );
if strcmp(method, 'ecc')
    % save eccentric points
    E0 = Ind;
    [tmp,a,E0] = unique(E0); E0 = reshape(E0,n,n);
    E0(M==0) = Inf;
    E0 = convert_distance_color(E0);
    for k=1:length(eccx)
        I = find( (X-eccx(k)).^2 + (Y-eccy(k)).^2 <= r^2 );
        for i=1:3
            A = E0(:,:,i); A(I) = i==1;
            E0(:,:,i) = A;
            A = E0(:,:,i); A(I) = i==1;
            E0(:,:,i) = A;
        end
    end
    imwrite( E0, [rep name '-shape-eccpoints.png'], 'png' );
    warning on;
end

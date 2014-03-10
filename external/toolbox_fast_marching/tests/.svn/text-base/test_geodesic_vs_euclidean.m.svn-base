% compare geodesic with euclidean distance inside a 2D shape


rep = 'images/';
name = 'bird';
name = 'giraffe';
name = 'camel';
name = 'chicken';

test = 'meangeodesic';
test = 'eccentricity';

% put s=-1 for Linf eccentricity
s = -1;

repsrc = 'data/';

rep = 'results/geodesic-vs-euclidean/';
if ~exist(rep)
    mkdir(rep);
end

n = 256;
M = rescale( load_image([repsrc name],n), 0,1 );
M = double(M>0.5);

% make sure pixels on the boundary are black
if M(1)==1
    M = 1-M;
end

% compare the geodesic ditance to the euclidean distance
clf;
imagesc(M); axis image; axis off;
[y,x] = ginput(1);
start_points = round([x y]');
W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;
options.constraint_map = [];
[D0,S,Q] = perform_fast_marching(W, start_points, options);
options.constraint_map = L;
[D1,S,Q] = perform_fast_marching(W, start_points, options);
D0(M==0) = Inf; D1(M==0) = Inf;

% compute a nice color image
D0a = convert_distance_color(D0);
D1a = convert_distance_color(D1);

% place a dot at the center
[Y,X] = meshgrid(1:n,1:n);
r = 0.02*n;
I = find( (X-x).^2 + (Y-y).^2 <= r^2 );
for i=1:3
    A = D0a(:,:,i); A(I) = i==1;
    D0a(:,:,i) = A;
    A = D1a(:,:,i); A(I) = i==1;
    D1a(:,:,i) = A;
end


% display
clf;
subplot(1,2,1)
imagesc(D0a); axis image; axis off; title('Euclidean');
subplot(1,2,2)
imagesc(D1a); axis image; axis off; title('Geodesic');

warning off;
imwrite(rescale(D0a), [rep name '-euclidean.png' ], 'png');
imwrite(rescale(D1a), [rep name '-geodesic.png' ], 'png');
warning on;

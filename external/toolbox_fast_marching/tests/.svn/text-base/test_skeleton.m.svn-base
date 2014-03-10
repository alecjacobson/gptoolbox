% test for skeleton and distance transform


rep  = 'data/';
name = 'mm';
name = 'cavern';
name = 'toto';

%% compute a binary shape
n = 200;
Ma = load_image([rep name],n-10);
Ma = sum(Ma,3);
M = zeros(n);
M(6:n-5,6:n-5) = Ma;
mask = 1-(M==M(1));

%% compute the skeleton
[skg,rad] = skeleton(mask);
% the higher the threshold, the smaller the skeleton
sk1 = skg>20;

sk2 = pick_curves(mask);

sk = {sk1,sk2};
nsk = length(sk);

clf;
for i=1:nsk
    %% compute euclidean distance to skeleton
    [Dsk1,Qsk] = eucdist2(logical(sk{i}));
    %% compute the geodesic distance, should not be exactly the same
    % positions along the skeleton
    [x,y] = ind2sub(size(mask),find(sk{i}));
    skpoints = [x(:)';y(:)'];
    CM = zeros(n) + Inf; CM(mask==0) = -Inf;
    options.constraint_map = CM; % constraint the propagation inside the shape
    [Dsk2,Z,Qsk] = perform_fast_marching(ones(n), skpoints, options);
    % display
    subplot(nsk,3,1+(i-1)*3);
    imagesc(mask+sk{i}); axis image; axis off;
    title('skeleton');
    subplot(nsk,3,2+(i-1)*3);
    imagesc(Dsk1.*mask); axis image; axis off;
    title('euclidean distance');
    subplot(nsk,3,3+(i-1)*3);
    imagesc(Dsk2.*mask); axis image; axis off;
    title('geodesic distance');
end
colormap jet(256);
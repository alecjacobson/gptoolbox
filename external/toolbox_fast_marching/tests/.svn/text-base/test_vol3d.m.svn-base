% volumetric data vizualization

% load brain1-128.mat
% load brain1-256.mat
load brain1-crop-256.mat


rep = 'results/vol3d/';
if not(exist(rep))
    mkdir(rep);
end

% save some slices
ilist = round(linspace(1,size(M,3),10)); ilist([1 end]) = [];
warning off;
for i=1:length(ilist)
    imwrite(rescale(M(:,:,ilist(i))), [rep 'brain-img-' num2str(i) '.png'], 'png');
end
warning on;


% crop
n = 100;
M = crop(M,n); % (end-n+1:end,end-n+1:end,end-n+1:end);
M = rescale(M);

clf;
h = vol3d('cdata',M,'texture','2D');
view(3); axis off;

center_list =  .55;% % .1:.1:.9;
nc = length(center_list);
for i=1:nc
    % set up alpha
    options.center = center_list(i);
    options.sigma = .08;
    a = compute_alpha_map('gaussian', options);
    colormap bone(256);
    vol3d(h);
%    saveas(gcf, [rep 'brain-vol3d-' num2str(i) '.png'], 'png');
end

n = size(M,1);
% perform fast marching
delta = 5;
clf;
imageplot(M(:,:,delta));
axis image; axis off;
title('Pick starting point');
start_point = round( ginput(1) );
start_point = [start_point(2); start_point(1); delta];


% compute the potential
c = M(start_point(1),start_point(2),start_point(3));
sigma = .04;
W = exp( -(M-c).^2 / (2*sigma^2) );
W = rescale(W,.001,1);
options.nb_iter_max = Inf;
fprintf('--> Perform Fast Marching ... ');
[D,S] = perform_fast_marching(W, start_point, options);
disp('done.');

% display the interface of the front
clf;
h = vol3d('cdata',rescale(D),'texture','2D');
view(3); axis off;
vol3d(h);

center_list =  .1; % [.1 .05 .02 .01 .005 .001];
options.sigma = 0.03;
nc = length(center_list);
for i=1:nc
    % set up alpha
    options.center = center_list(i);
    options.cname = 'jet';
    [a,c] = compute_alpha_map('gaussian', options);
    colormap(c);
    vol3d(h);
end

% find an end point
d = D(:,:,n-delta);
[tmp,I] = min(d(:));
[x,y] = ind2sub(size(d),I);
end_point = [x;y;n-delta];

path = compute_discrete_geodesic(D,end_point);

% draw the path
Dend = D(end_point(1),end_point(2),end_point(3));
D1 = double( D<=Dend );
clf;
plot_fast_marching_3d(M,D1,path,start_point,end_point);
view(3);  camlight;


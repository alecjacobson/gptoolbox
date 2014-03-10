% shows an animation of a propagation
%
%   Copyright (c) 2005 Gabriel Peyre

path(path, 'toolbox/');
n = 100;
name = 'road2';
name = 'peaks';
name = 'mountain';
name = 'constant';
[M,W] = load_potential_map(name, n);
W = rescale(W,0.01,1);


[start_points,end_points] = pick_start_end_point(W);

k = 500*n/400;
alpha = 1.7;
klist = round( k*alpha.^(0:log2(n^2/k)/log2(alpha)) ); klist(end+1) = n^2;
klist = Inf;

save_images = 1;
rep = 'results/propagation-2d/';
if not(exist(rep))
    mkdir(rep);
end    

clf;
colormap gray(256);
for i=1:length(klist)
    progressbar(i,length(klist));
    options.nb_iter_max = klist(i);
    [D,S] = perform_fast_marching(W, start_points, options);
    % plot_fast_marching_2d(W, S);
    % compute a nice color image
    c = jet(256);
    U = D; U(U==Inf) = 0;
    I = floor(255*rescale(U))+1;
    v = c(I(:),:); v = reshape(v, [n n 3]);
    A = repmat(rescale(M), [1 1 3]);
    B = repmat(D, [1 1 3]);
    v(B==Inf) = A(B==Inf);
    % display
    subplot(3, ceil(length(klist)/3), i);
    imagesc(v); axis image; axis off;
    if save_images
        warning off;
        imwrite(rescale(v), [rep name '-propagation-' num2string_fixeddigit(i, 3) '.png' ], 'png');
        warning on;
    end
end
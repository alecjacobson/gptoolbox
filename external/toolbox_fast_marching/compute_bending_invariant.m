function vertex1 = compute_bending_invariant(vertex,faces,options)

% compute_bending_invariant - compute bending invariants
%
%   vertex1 = compute_bending_invariant(vertex,faces,options);
%
%   Use out-of-sample interpolation to speed up the computation.
%   The number of landmarks used is set by options.nlandmarks.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
nlandmarks = getoptions(options, 'nlandmarks', 200);
nverts = max(size(vertex));
nlandmarks = min(nlandmarks,nverts);

% use farthest point sampling to compute the landmarks
Dland = [];
for i=1:nlandmarks
    if i==1
        landmarks = floor(rand*nverts)+1;
    else
        [tmp,landmarks(end+1)] = max( min(Dland, [], 2) );
    end
    progressbar(i,nlandmarks);
    [d,S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(i));
    Dland(:,end+1) = d(:);
end
fprintf('\n');
% square distances
Dland = Dland.^2;
% perform isomap on the reduced set of points
D1 = Dland(landmarks,:); % reduced pairwise distances
D1 = (D1+D1')/2; % force symmetry
J = eye(nlandmarks) - ones(nlandmarks)/nlandmarks; % centering matrix
K = -1/2 * J*D1*J; % inner product matrix
% compute the rank-3 approximation of the inner product to compute embedding
opt.disp = 0;
[xy, val] = eigs(K, 3, 'LM', opt);
xy = xy .* repmat(1./sqrt(diag(val))', [nlandmarks 1]);
% interpolation on the full set of points
% extend the embeding using geodesic interpolation
vertex1 = zeros(nverts,3);
deltan = mean(Dland,1);
for x=1:nverts
    deltax = Dland(x,:);
    vertex1(x,:) = 1/2 * ( xy' * ( deltan-deltax )' )';
end

return;

clf;
hold on;
plot_mesh(vertex, faces, options);
shading interp; camlight;
h = plot3(v(1,:), v(2,:), v(3,:), '.');
set(h, 'MarkerSize', 20);
hold off;
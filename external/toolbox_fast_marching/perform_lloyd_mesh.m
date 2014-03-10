function start_points = perform_lloyd_mesh(vertex,faces, start_points, options)

% perform_lloyd_mesh - perform lloyd relaxation to sample point on a mesh
%
%   start_points = perform_lloyd_mesh(vertex,faces, start_points, options);
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end
if size(faces,1)>size(faces,2)
    faces = faces';
end

niter_lloyd = getoptions(options, 'niter_lloyd', 1);
if niter_lloyd>1
    for i=1:niter_lloyd
        start_points = perform_lloyd_mesh(vertex,faces, start_points, options)
        options.lambda = []; % enfore recomputing
    end
    return;
end

lambda = getoptions(options, 'lambda', []);
edges_id = getoptions(options, 'edges_id', []);
Q = getoptions(options, 'Q', []);
if isempty(lambda) || isempty(edges_id) || isempty(Q)
    % update voronoi
    [Q,DQ, ve, edges_id, lambda] = compute_voronoi_mesh(vertex,faces, start_points, options);
end

% compute distances for start
ne = length(lambda);
n = size(vertex,2);
% compute edge length
d = diff( reshape(vertex(:,edges_id'), [3 ne 2]), 1, 3 );
d = sqrt( sum(d.^2,1) );
d1 = (1-lambda).*d;
d2 = lambda.*d;
% compute initial values for FM
values = zeros(n,1); cnt = zeros(n,1);
for i=1:ne
    values(edges_id(1,i)) = values(edges_id(1,i)) + d1(i);
    cnt(edges_id(1,i)) = cnt(edges_id(1,i)) + 1;
    values(edges_id(2,i)) = values(edges_id(2,i)) + d1(i);
    cnt(edges_id(2,i)) = cnt(edges_id(2,i)) + 1;
end
values = values./cnt;
% perform FM from edge points
pts = unique( edges_id(:) );
options.values = values(pts);
options.values = [];
D = perform_fast_marching_mesh(vertex, faces, pts, options);

% seed new locations
nstart = length(start_points);
for i=1:nstart
    I = find(Q(:,1)==i);
    [v,k] = max(D(I));
    start_points(i) = I(k);
end
start_points = start_points(:);
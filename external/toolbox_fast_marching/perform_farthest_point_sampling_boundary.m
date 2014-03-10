function points = perform_farthest_point_sampling_boundary( W, points, nbr_pts, method )

% perform_farthest_point_sampling - samples points using farthest seeding strategy
%
% points = perform_farthest_point_sampling( W, points, nbr_points );
%
%   points can be []
%   
%   Copyright (c) 2005 Gabriel Peyré

if nbr_pts>1
    for i=1:nbr_pts
        points = perform_farthest_point_sampling_boundary( W, points, 1 );
    end
    return;
end

if isempty(points)
    points = [1;1];
    return;
end

if size(points,1)~=2
    points = points';
end
n = size(W,1);

if nargin<4
    method = 'boundary';
end

if strcmp(method, 'boundary')
    t = 1:n;
    it = 1:n;
    a = 0*t+1;
    b = 0*t+n;
    trial = [ [t;a] [b;t] [it;b] [a;it] ];
else
    x = round( linspace(1,n,64) );
    [Y,X] = meshgrid(1:n,1:n);
    trial = [X(:)'; Y(:)'];
end

D = compute_distance_to_points(points,trial);

d = min(D,[],2);

[tmp,I] = max(d);

points = [points, trial(:,I)];
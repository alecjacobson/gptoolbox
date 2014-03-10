function D = compute_distance_to_points(X,seeds)

% compute_distance_to_points - compute euclidean distance to a set of points.
%
%   D = compute_distance_to_points(X,seeds)
%
%   'X' is a [d,n] matrix, X(:,i) is the ith point living in R^d.
%   'seeds' is a [d,k] matrix.
%   D(i,j) = |X(:,j)-seeds(:,i)|^2.
%
%   Copyright (c) 2004 Gabriel Peyré

nbCluster = size(seeds,2);
n = size(X,2);
D = zeros(nbCluster,n);
d = size(X,1);

for k=1:nbCluster
    % distance to seed
    D(k,:) = sum( (X - repmat(seeds(:,k),1,n)).^2 );
end
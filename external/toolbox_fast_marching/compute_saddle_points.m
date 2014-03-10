function [vertex,faces] = compute_saddle_points(Q,D,mask)

% compute_saddle_points - compute saddle points of a Voronoi segmentation
%
%   [vertex,faces] = compute_saddle_points(Q,D[,mask]);
%
%   Q is a voronoi index map.
%   D is the distance function associated to the voronoi map.
%
%   The vertex of the saddle points are first the double points (meeting
%   points of the voronoi diagram along the boundary of the domain) and
%   then the triple points (meeting points of 3 cells).
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin==3 && not(isempty(mask))
    Q(mask==0) = -1;
end

Q1 = zeros(size(Q)+2)-1;
Q1(2:end-1,2:end-1) = Q;
V = [];
v = Q1(1:end-1,1:end-1); V = [V v(:)];
v = Q1(2:end,1:end-1); V = [V v(:)];
v = Q1(1:end-1,2:end); V = [V v(:)];
v = Q1(2:end,2:end); V = [V v(:)];
V = sort(V,2);
d = (V(:,1)~=V(:,2)) + (V(:,2)~=V(:,3)) + (V(:,3)~=V(:,4));
V = V';

I = find(d>=2);

[vx,vy] = ind2sub(size(Q)+1, I);
vx = clamp(vx,1,size(Q,1));
vy = clamp(vy,1,size(Q,1));
J = vx+(vy-1)*size(Q,1);

% sort according to distance
[tmp,s] = sort(D(J), 1, 'descend');
I = I(s);

[vx,vy] = ind2sub(size(Q)+1, I);
vx = clamp(vx,1,size(Q,1));
vy = clamp(vy,1,size(Q,1));
vertex = cat(1,vx',vy');

V = sort(V, 1, 'descend');
faces = V(1:3, I);

if isempty(vertex)
    % add farthest point
    [tmp,I] = max( D(:) );
    [vx,vy] = ind2sub(size(D), I(1));
    vertex = [vx;vy];
    faces = [-1 -1 -1]';
end

I = find( faces(1,:)<0 | faces(2,:)<0 | faces(3,:)<0 );
J = find( faces(1,:)>0 & faces(2,:)>0 & faces(3,:)>0 );
faces = cat(2, faces(:,I), faces(:,J) );

function h = plot_edges(edges, vertex, color)

% plot_edges - plot a list of edges
%
%   h = plot_edges(edges, vertex, color);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    str = 'b';
end

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end

x = [ vertex(1,edges(1,:)); vertex(1,edges(2,:)) ];
y = [ vertex(2,edges(1,:)); vertex(2,edges(2,:)) ];
if size(vertex,1)==2
    h = line(x,y, 'color', color);
elseif size(vertex,1)==3
    z = [ vertex(3,edges(1,:)); vertex(3,edges(2,:)) ];
    h = line(x,y,z, 'color', color);    
else
    error('Works only for 2D and 3D plots');    
end
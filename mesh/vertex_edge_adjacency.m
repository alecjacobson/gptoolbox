function VE = vertex_edge_adjacency(E)
%VERTEX_EDGE_ADJACENCY Build a vertex-edge adjacency matrix for a triangle
% mesh.
%
% Input:
%   E   #E x 2   list of edge indices
% Output:
%   VE  #E x #V  sparse matrix, VE(i,j) is 1 iff vertex j is in edge i
%
% Example:
%   [EF,EI,uE,EMAP] = edge_flaps(F);
%   VE = vertex_edge_adjacency(uE);
i = (1:size(E,1))';
j = E;
VE = sparse([i i],j,1);
end
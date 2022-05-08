function [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH)
% BAD_QUAD_MESH_FROM_QUADTREE
% From a proper quadtree, builds a connected but degenerate quad mesh
% containing only the leaf nodes, mostly for visualization purposes.
%
% [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH)
%
% Inputs:
%   C #nodes by 3 matrix of cell centers
%   W #nodes vector of cell widths (**not** half widths)
%   CH #nodes by 4 matrix of child indeces (-1 if leaf node)
%
% Outputs:
%   V #V by 3 matrix of vertex positions
%   Q #Q by 4 matrix of quad indeces into V
%
% Example:
%
% th = 2*pi*rand(200,1);
% P = 0.5*[cos(th),sin(th)];
% P = [P;[-1,-1];[1,1]];
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'MaxDepth',8,'Graded',true);
% [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH);
% 
% clf
% hold off
% hold on
% tsurf(Q,V,'FaceColor','b','EdgeColor','k',falpha(0.2,1))
% axis equal
%
% See also: initialize_quadtree.m

is_child = (CH(:,1)==-1);
W = W(is_child,:);
C = C(is_child,:);
V = [C + 0.5*W*[-1,-1];C + 0.5*W*[-1,1];C + 0.5*W*[1,1];C + 0.5*W*[1,-1]];
Q = [1:size(W,1)]';
Q = [Q,Q+size(W,1),Q+2*size(W,1),Q+3*size(W,1)];
[V,SVI,SVJ] = remove_duplicate_vertices(V,min(W)/100);
% remap faces
Q = SVJ(Q);
% I hate this part
if size(Q,2)==1
    Q = Q';
end


end
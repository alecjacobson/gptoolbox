% DECIMATE_CGAL Decimate a mesh (V,F) using CGAL's Polyhedron edge_collapse
% routine. (V,F) must already be a manifold mesh
%
% [W,G] = decimate_cgal(V,F,ratio)
%
% Inputs:
%   V  #V by 3 list of input mesh vertex positions
%   F  #F by 3 list of triangles indices into V
%   ratio  scalar between 0 and 1 exclusively specifying ratio of output faces
%   to input faces: ratio ~= #G / #F
% Outputs:
%   W  #W by 3 list of output mesh vertex positions
%   G  #G by 3 list of triangles indices into W
%
% See also: reducepatch
%

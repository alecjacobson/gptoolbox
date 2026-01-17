% RAY_MESH_INTERSECT_ALL Find all hits for each ray.
%
% [I, J, T, lambda] = ray_mesh_intersect_all(source, dir, V, F)
%
% Input:
%    source #rays by 3 list of 3D vector ray origins
%    dir #rays by 3 list of 3D vector ray directions
%    V  #V by 3 list of vertex positions
%    F  #F by 3 list of triangle indices
% Output:
%    I  #hits list of indices into 1:#rays
%    J  #hits list of indices into 1:#F
%    T  #hits list of parametric distances along the rays
%    lambda   #hits by 3 list of barycentric coordinates
%


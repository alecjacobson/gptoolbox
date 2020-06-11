% RAY_MESH_INTERSECT  Find first hit (if it exists) for each ray.
%
% [flag, t, lambda] = ray_mesh_intersect(src, dir, V, F);
%
% Input:
%    src #rays by 3 list of 3D vector ray origins
%    dir #rays by 3 list of 3D vector ray directions
%    V  #V by 3 list of vertex positions
%    F  #F by 3 list of triangle indices
% Output:
%    id  #rays list of indices into F (0 if no hit)
%    t  #rays list of distances from the ray origin (inf if no hit)
%    lambda  #rays by 3 list of barycentric coordinate of hit
%

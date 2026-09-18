% REFINE_TRIANGULATION Refine an existing constrained Delaunay triangulation
% via triangle via libigl.
%
% [TV,TF] = refine_triangulation(V,E,F);
% [TV,TF] = refine_triangulation(V,E,F,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   V  #V by 2 list of input vertex positions
%   E  #E by 2 list of vertex ids forming segments (may be empty)
%   F  #F by 3 list of input triangles as indices into rows of V
%   Optional:
%     'Flags'  followed by additional flags to pass to triangle (the 'r'
%       reconstruct flag and 'z' zero-indexing are always added) {''}
% Outputs:
%   TV  #TV by 2 list of output vertex positions (V should appear as first
%     rows)
%   TF  #TF by 3 list of output triangles as indices into rows of TV
%
% See also: triangulate

% TRIANGULATE Constrained Delaunay triangulation via triangle via libigl.
%
% [TV,TF,TVM,TEM] = triangulate(V,E);
% [TV,TF,TVM,TEM] = triangulate(V,E,'ParameterName',ParameterValue, ...)
%
% The overhead of calling triangle via the command line can be significant. For
% example, to trivially triangulate three points. Calling `triangulate` is
% 15000x faster than `triangle` on macos 10.15 and Matlab 2020b.
%
% Inputs:
%   V  #V by 2 list of points to triangulate
%   E  #E by 2 list of edges as indices into rows of V
%   Optional:
%     'EdgeMarkers'  followed by #E list of edge markers
%     'Flags'  followed by flags
%     'Holes' followed by #H by 2 list of hole positions
%     'VertexMarkers'  followed by #V list of vertex markers
% Outputs:
%   TV  #TV by 2 list of output points (V should appear as first rows)
%   TF  #TF by 3 list of triangles as indices into rows of TV
%   TVM  #TV list of output vertex markers
%   TE  #TE list of output segments (inputs and boundary)
%   TEM  #TE list of output edge markers. Implicit boundaries (e.g., from "-c")
%     will be marked as max(EM)+1

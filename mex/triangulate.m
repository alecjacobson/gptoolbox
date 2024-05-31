% TRIANGULATE Constrained Delaunay triangulation via triangle via libigl.
%
% [TV,TF,TVM,TE,TEM] = triangulate(V,E);
% [TV,TF,TVM,TE,TEM] = triangulate(V,E,'ParameterName',ParameterValue, ...)
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
%
% Known issues:
%   igl doesn't expose the 'r' Refine flag

% Example:
%   % Work around for meshing with holes without over refining near small holes:
%   [V,F] = triangulate(OV,OE);
%   BC = barycenter(V,F);
%   H = BC(abs(winding_number(V,E,BC))<0.5,:);
%   [V,F] = triangulate(V,E,'Holes',H,'Flags','-q');
%
%   % compare to:
%   [V2,F2] = triangulate(OV,OE,'Flags','-q33');
%   F2 = F2(abs(winding_number(OV,OE,barycenter(V2,F2)))>0.5,:);
%   [V2,~,~,F2] = remove_unreferenced(V2,F2);
%   tsurf(F,V)
%   hold on;
%   tsurf(F2,V2+[max(V(:,1))-min(V(:,1)) 0]);
%   hold off;
%

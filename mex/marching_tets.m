% [SV,SF,J,BC] = marching_tets(TV,TT,S)
% [SV,SF,J,BC] = marching_tets(TV,TT,S,'ParameterName',ParameterValue, ...)
%
% Compute the isosurface of a scalar field defined on the vertices of a
% tetrahedral mesh using the marching tetrahedra algorithm.
%
% Inputs:
%   TV  #TV by 3 list of tetrahedral mesh vertex positions
%   TT  #TT by 4 list of tetrahedra indices into rows of TV
%   S  #TV by 1 list of scalar values at each vertex of TV
%   Optional:
%     'IsoValue' followed by the isovalue of the level set to compute {0}
% Outputs:
%   SV  #SV by 3 list of output level-surface mesh vertex positions
%   SF  #SF by 3 list of output level-surface mesh triangle indices into
%     rows of SV
%   J  #SF list of indices into rows of TT revealing which tet each face
%     comes from
%   BC  #SV by #TV sparse matrix of barycentric coordinates so that
%     SV = BC*TV
%
% See also: isolines, isolines_intrinsic

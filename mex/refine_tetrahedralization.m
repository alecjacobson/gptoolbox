% REFINE_TETRAHEDRALIZATION Refine an existing tetrahedralization via tetgen
% via libigl.
%
% [RV,RT,RF] = refine_tetrahedralization(TV,TT);
% [RV,RT,RF] = refine_tetrahedralization(TV,TT,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   TV  #TV by 3 list of input vertex positions
%   TT  #TT by 4 list of input tetrahedra as indices into rows of TV
%   Optional:
%     'Flags'  followed by additional flags to pass to tetgen (the 'r'
%       reconstruct flag and 'Q' quiet flag are always added) {''}
%     'Boundary' followed by #TF by 3 list of known boundary triangles as
%       indices into rows of TV, used to preserve the input boundary
%       {computed automatically via boundary_facets(TT)}
% Outputs:
%   RV  #RV by 3 list of output vertex positions (TV should appear as first
%     rows)
%   RT  #RT by 4 list of output tetrahedra as indices into rows of RV
%   RF  #RF by 3 list of output boundary triangles as indices into rows of
%     RV
%
% See also: tetrahedralize, refine_triangulation

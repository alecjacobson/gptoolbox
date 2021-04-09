% TETRAHEDRALIZE
%
% [TV,TT] = tetgen(SV,SF)
% [TV,TT,TF,TR,TN,PT,FT] = tetgen(SV,SF,'Flags',flags)
%
% Inputs:
%   SV  #SV by 3 list of input vertex positions
%   SF  #SF by 3 list of face indices into rows of SV
%   Optional:
%     'Flags'  followed by TetGen flags to use {'-q2'}
% Outputs:
%   TV  #TV by 3 list of output vertex positions (SV should appear as first
%     rows)
%   TT  #TT by 4 list of output tetrahedron indices into rows of TV
%   TF  #TF by 3 list of output boundary triangles indices into rows of TV
%   TR  #TT list of region ID for each tetrahedron      
%   TN  #TT by 4 list of indices neighbors for each tetrahedron
%   PT  #TV list of incident tetrahedron for a vertex
%   FT  #TF by 2 list of tetrahedrons sharing a triface      
%   
% 
%

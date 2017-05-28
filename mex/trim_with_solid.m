% TRIM_WITH_SOLID Given an arbitrary mesh (VA,FA) and the boundary mesh (VB,FB)
% of a solid (as defined in [Zhou et al. 2016]), Resolve intersections between
% A and B subdividing faces of A so that intersections with B exists only along
% edges and vertices (and coplanar faces). Then determine whether each of these
% faces is inside or outside of B. This can be used to extract the part of A
% inside or outside of B.
%
% Inputs:
%   VA  #VA by 3 list of mesh vertex positions of A
%   FA  #FA by 3 list of mesh triangle indices into VA
%   VB  #VB by 3 list of mesh vertex positions of B
%   FB  #FB by 3 list of mesh triangle indices into VB
% Outputs:
%   V  #V by 3 list of mesh vertex positions of output
%   F  #F by 3 list of mesh triangle indices into V
%   D  #F list of bools whether face is inside B
%   J  #F list of indices into FA revealing birth parent
%

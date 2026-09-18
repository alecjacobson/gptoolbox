% [iB,iFI,iE,I] = isolines_intrinsic(F,S,vals)
%
% Compute isolines of a scalar field on a triangle mesh intrinsically (i.e.,
% without reference to vertex positions).
%
% Inputs:
%   F  #F by 3 list of mesh triangle indices into some V
%   S  #S by 1 list of per-vertex scalar values
%   vals  #vals by 1 list of values to compute isolines for
% Outputs:
%   iB  #iB by 3 list of barycentric coordinates so that
%     iV(i,:) = iB(i,1)*V(F(iFI(i),1),:) + ...
%               iB(i,2)*V(F(iFI(i),2),:) + ...
%               iB(i,3)*V(F(iFI(i),3),:)
%   iFI  #iB list of triangle indices for each row of iB (all points will
%     either lie on an edge or vertex: an arbitrary incident face will be
%     given)
%   iE  #iE by 2 list of edge indices into rows of iB
%   I  #iE by 1 list of indices into vals indicating which value each
%     segment belongs to
%
% See also: isolines

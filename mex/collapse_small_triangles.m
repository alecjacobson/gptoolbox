% Given a triangle mesh (V,F) compute a new mesh (VV,FF) which contains the
% original faces and vertices of (V,F) except any small triangles have been
% removed via collapse.
%
% We are *not* following the rules in "Mesh Optimization" [Hoppe et al]
% Section 4.2. But for our purposes we don't care about this criteria.
%
% [FF] = collapse_small_triangles(V,F,eps)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into V
%   eps  epsilon for smallest allowed area treated as fraction of squared bounding box
%     diagonal
% Outputs:
%   FF  #FF by 3 list of triangle indices into V
%

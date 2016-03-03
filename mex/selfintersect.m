% SELFINTERSECT Given a triangle mesh (V,F) compute a new mesh (VV,FF) which is the same as
% (V,F) except that any self-intersecting triangles in (V,F) have been
% subdivided (new vertices and face created) so that the self-intersection
% contour lies exactly on edges in (VV,FF). New vertices will appear in
% original faces or on original edges. New vertices on edges are "merged" only
% across original faces sharing that edge. This means that if the input
% triangle mesh is a closed manifold the output will be too.
% 
% [VV,FF,IF] = selfintersect(V,F,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into V
%   Optional:
%     'DetectOnly'  followed by bool. Whether to only detect intersecting pairs
%       (sets IF but not VV or FF) {false}
%     'FirstOnly'  followed by bool whether to only detect the first
%       intersection. {false}
%     'StitchAll'  followed by whether to stitch all vertices in the output, if
%       true then IM will be 1:size(VV,1) {false}
% Outputs:
%   VV  #VV by 3 list of vertex positions
%   FF  #FF by 3 list of triangle indices into V
%   IF  #intersecting face pairs by 2  list of intersecting face pairs,
%     indexing F
%   J   #FF list of indices into F of birth parents
%   IM  #VV list of indices into VV of unique vertices
%
% Example:
%   [SV,SF,~,~,IM] = selfintersect(V,F);
%   FF = IM(SF);
%   [U,IM] = remove_unreferenced(SV,FF);
%   G = IM(FF);
%
%   % Self-Intersect (V,F)+(U,G) and separate
%   VU = [V;U];
%   FG = [F;size(V,1)+G];
%   [SVU,SFG,~,J] = selfintersect(VU,FG);
%   SF = SFG(J<=size(F,1),:);
%   SG = SFG(J>size(F,1),:);
%   [SV,IM] = remove_unreferenced(SVU,SF);
%   SF = IM(SF);
%   [SU,IM] = remove_unreferenced(SVU,SG);
%   SG = IM(SG);
%
%

% See selfintersect.h, selfintersect.cpp for mex implementation

% INTERSECT_OTHER Given a triangle mesh (V,F) and another mesh (U,G) find all
% pairs of intersecting faces. Note that self-intersections are ignored.
% 
% [IF] = intersect_other(V,F,U,G,'ParameterName',ParameterValue, ...)
% [IF,VVAB,FFAB,JAB,IMAB] = intersect_other(V,F,U,G,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into V
%   U  #U by 3 list of vertex positions
%   G  #G by 3 list of triangle indices into U
%   Optional:
%     'FirstOnly'  followed by bool whether to only detect the first
%       intersection.
% Outputs:
%   IF  #intersecting face pairs by 2 list of intersecting face pairs,
%     indexing F and G
%   VVAB  #VVAB by 3 list of vertex positions
%   FFAB  #FFAB by 3 list of triangle indices into VVA
%   JAB  #FFAB list of indices into [FA;FB] denoting birth triangle
%   IMAB  #VVAB list of indices stitching duplicates (resulting from
%     mesh intersections) together
%
% See also: selfintersect

% See intersect_other.h, intersect_other.cpp for mex implementation


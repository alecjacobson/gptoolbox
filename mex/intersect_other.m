% INTERSECT_OTHER Given a triangle mesh (V,F) and another mesh (U,G) find all
% pairs of intersecting faces. Note that self-intersections are ignored.
% 
% [IF] = intersect_other(V,F,U,G,'ParameterName',ParameterValue, ...)
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
%
% See also: selfintersect

% See intersect_other.h, intersect_other.cpp for mex implementation


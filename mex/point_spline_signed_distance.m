% [signedD,I,S,K,T] = point_spline_signed_distance(Q,P,C)
% [signedD,I,S,K,T,B1,B2,leaf] = point_spline_signed_distance(Q,P,C,'Acceleration',{B1,B2,leaf});
%
% Inputs:
%   Q   #Q by dim list of query points
%   P   #P by dim list of control points
%   C   #C by 4 list of indices into P of cubic splines
%   Optional:
%     'Acceleration' followed by cell containing B1, B2, leaf
% Outputs:
%   signedD #Q list of signed distances distances
%   I  #Q list of indices into C
%   S  #Q list of parameter value on corresponding cubic spline
%   K  #Q by dim list of closest points on corresponding cubic spline
%   T  #Q by dim list of tangent vectors on corresponding cubic spline
%   B1  #B by dim list of min corners of eytzinger AABB tree boxes
%   B2  #B by dim list of max corners of eytzinger AABB tree boxes
%   leaf  #B list of indices in to C/flags
%   

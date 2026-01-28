% SPLINE_WINDING_NUMBER Compute generalized winding numbers for a query
% points O with respect to a given cubic Bezier spline (P,C).
%
% W = spline_winding_number(P,C,O)
% [W,B1,B2,leaf] = spline_winding_number(P,C,O,'Acceleration',{B1,B2,leaf})
%
% Inputs:
%   P  #P by 2 list of control point locations
%   C  #C by 4 list of indices into P of cubic Bezier curves
%   O  #O by 2 list of query points
%   Optional:
%     'Acceleration' followed by cell containing B1, B2, leaf
% Outputs:
%   W  #O list of winding number values
%   B1  #B by dim list of min corners of eytzinger AABB tree boxes
%   B2  #B by dim list of max corners of eytzinger AABB tree boxes
%   leaf  #B list of indices in to C/flags
%

function [ A ] = internalangles_intrinsic(l)
% INTERNALANGLES_INTRINSIC
%
% A = internalangles_intrinsic(l)
%
% Compute internal angles per face (in degrees) using edge lengths
%
% Inputs:
%  l  #F x 3  matrix of edge lengths
% Output:
%  A  #F x 3 list of triples of triangle angles
%

s12 = l(:,3);
s13 = l(:,2);
s23 = l(:,1);

a12 = acos((s13.^2 + s23.^2 - s12.^2)./(2.*s13.*s23));
a13 = acos((s12.^2 + s23.^2 - s13.^2)./(2.*s12.*s23));
a23 = acos((s12.^2 + s13.^2 - s23.^2)./(2.*s12.*s13));

A = [a23 a13 a12];

end



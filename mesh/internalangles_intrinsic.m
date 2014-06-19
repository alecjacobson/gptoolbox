function [ A ] = internalangles_intrinsic(l)
% INTERNALANGLES_INTRINSIC Compute internal angles per face (in degrees) using
% edge lengths
%
% A = internalangles_intrinsic(l)
%
% Inputs:
%  l  #F x 3  matrix of edge lengths
% Output:
%  A  #F x 3 list of triples of triangle angles
%

s23 = l(:,1);
s31 = l(:,2);
s12 = l(:,3);

a23 = acos((s12.^2 + s31.^2 - s23.^2)./(2.*s12.*s31));
a31 = acos((s23.^2 + s12.^2 - s31.^2)./(2.*s23.*s12));
a12 = acos((s31.^2 + s23.^2 - s12.^2)./(2.*s31.*s23));

A = [a23 a31 a12];

end



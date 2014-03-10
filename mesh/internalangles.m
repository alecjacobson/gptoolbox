function [ A ] = internalangles( V, F)
% INTERNALANGLES Compute internal angles per face (in degrees)
%
% A = internalangles(V,F)
%
% Inputs:
%  V  #V x 3 matrix of vertex coordinates
%  F  #F x 3  matrix of indices of triangle corners
% Output:
%  A  #F x 3 list of triples of triangle angles
%
i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);

s12 = normrow(V(i2,:) - V(i1,:));
s13 = normrow(V(i3,:) - V(i1,:));
s23 = normrow(V(i3,:) - V(i2,:));

l = [s23 s13 s12];

A = internalangles_intrinsic(l);

end


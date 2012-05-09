function [ dblA ] = doublearea( V, F )
  % DOUBLEAREA Compute the double area of the triangles of a mesh
  %
  % [ dblA ] = doublearea( V, F )
  %
  %
  % Inputs:
  %  V #V x 3 matrix of vertex coordinates
  %  F #F x 3  matrix of indices of triangle corners
  % Outputs:
  %  dblA   #F list of twice the area of each corresponding face
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), and Daniele Panozzo
  %
  
  % triangles
  % edge lengths numbered same as opposite vertices
  l = [ ...
    sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
    sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
    sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
    ];
  l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
  % semiperimeters
  s = (l1 + l2 + l3)*0.5;
  % Heron's formula for area
  dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));

end

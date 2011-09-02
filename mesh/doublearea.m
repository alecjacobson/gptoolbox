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
  
  i1 = F(:,1);
  i2 = F(:,2);
  i3 = F(:,3);
  
  % append zeros so cross product will work correctly for 2D vertices
  if(size(V,2) == 2)
    V = [V zeros(size(V,1),1)];
  end

  T = cross(V(i2,:)-V(i1,:),V(i3,:)-V(i1,:),2);
  
  dblA = sqrt(sum(T.^2,2)); % Face areas
  
end


function [ N ] = normals( V, F )
  % NORMALS 
  %
  % N = normals(V,F)
  %
  % Compute *unnormalized* normals per face
  %
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x 3  matrix of indices of triangle corners
  % Output:
  %  N  #F x 3 list of face normals
  %
  
  p1 = V(F(:,1),:);
  p2 = V(F(:,2),:);
  p3 = V(F(:,3),:);
  
  N = cross(p2 - p1, p3 - p1);

end


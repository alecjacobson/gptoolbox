function [ B ] = barycenter(V, F)
  % BARYCENTER
  %
  % B = barycenter(V,F)
  %
  % Compute the barycenter of every triangle
  %
  % Inputs:
  %   V #V x 3 matrix of vertex coordinates
  %   F #F x 3  matrix of indices of triangle corners
  % Output:
  %   B a #F x 3 matrix of 3d vertices

  i1 = F(:,1);
  i2 = F(:,2);
  i3 = F(:,3);
  
  B = (V(i1,:) + V(i2,:) + V(i3,:)) / 3.0;

end


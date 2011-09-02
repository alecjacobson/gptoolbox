function [ b ] = avgedge(V,F)
  % AVGEDGE Compute the average of every edge in the mesh
  % 
  % [ b ] = avgedge(V,F)
  %
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x 3  matrix of indices of triangle corners
  % Outputs:
  %  b average edge length
  %
  % Note: boundary edges are weighted half as much as internal edges
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), and Daniele Panozzo
  %

  
  i1 = F(:,1);
  i2 = F(:,2);
  i3 = F(:,3);
  
  B = (                                ...
          normrow(V(i1,:) - V(i2,:)) + ...
          normrow(V(i2,:) - V(i3,:)) + ...
          normrow(V(i1,:) - V(i3,:))   ...
      ) / 3.0;
  
  b = mean(B);

end


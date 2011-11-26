function [ b ] = avgedge(V,F)
  % AVGEDGE Compute the average of every edge in the mesh
  % 
  % [ b ] = avgedge(V,F)
  %
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x #simplex size  list of simplex indices
  % Outputs:
  %  b average edge length
  %
  % Note: boundary edges are weighted half as much as internal edges
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), and Daniele Panozzo
  %

  
  % works on anything edges.m can handle
  E = edges(F);
  B = normrow(V(E(:,1),:)-V(E(:,2),:));
  
  b = mean(B);

end


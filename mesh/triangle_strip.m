function [V,F] = strip(P,E)
  % STRIP build a strip of triangles given a curve
  %
  % [V,F] = strip(P,E)
  % 
  % Inputs:
  %   P  #points by 2 list of 2d points
  %   E  #edges by 2 list of edge indices into P
  % Outputs:
  %   V  #points*2 by 3 list of 3d points
  %   F  #edges*2 by 3 list of triangle indices
  %
  n = size(P,1);
  V = [P -ones(n,1);P ones(n,1)];
  F = [E n+E(:,1);n+fliplr(E) E(:,2)];
end

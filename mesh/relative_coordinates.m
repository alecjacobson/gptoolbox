function [X,Y] = relative_coordinates(V,F)
  % RELATIVE_COORDINATES computes the relative coordinates of each triangle
  % corner. That is, each triangle coordinate can be written in the basis of
  % the opposite edge and that edge rotated 90 degrees. Thus each corner can be
  % written in terms of the other two corners.
  %
  % [X,Y] = relative_coordinates(V,F) 
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   X  #F by 3 list of coordinates corresponding to each corners opposite
  %     edge
  %   Y  #F by 3 list of coordinates corresponding to each corners opposite
  %     edge, rotated 90 degrees
  %

  % Indices of each triangle vertex, I, and its corresponding two neighbors, J
  % and K
  I = [F(:,1);F(:,2);F(:,3)];
  J = [F(:,2);F(:,3);F(:,1)];
  K = [F(:,3);F(:,1);F(:,2)];
  % rename for rest positions of I,J,K for "convenience"
  % x-coordinates
  v0x = V(J,1);
  v1x = V(K,1);
  v2x = V(I,1);
  % y-coordinates
  v0y = V(J,2);
  v1y = V(K,2);
  v2y = V(I,2);
  % Compute, X and Y so that we can represent V(I,:) in terms of V(J,:) and
  % V(K,:), namely:
  % 
  % V(I,:) = 
  %   (V(J,:) + [X X].* (V(K,:) - V(J,:)) + [Y -Y].* fliplr(V(K,:)-V(J,:)))
  %
  X = ...
    (v0x.^2 - v1x.*v0x - v0x.*v2x + v1y.*v2y - v1y.*v0y - v0y.*v2y + v2x.*v1x + v0y.^2)./ ...
    (v1x.^2 - 2.*v1x.*v0x + v0x.^2 + v1y.^2 - 2.*v1y.*v0y + v0y.^2);
  Y = ...
    -(v1x.*v2y - v1x.*v0y - v1y.*v2x + v1y.*v0x - v0x.*v2y + v0y.*v2x)./ ...
    (v1x.^2 - 2.*v1x.*v0x + v0x.^2 + v1y.^2 - 2.*v1y.*v0y + v0y.^2);

  % reshape to match F
  X = reshape(X,size(F,1),3);
  Y = reshape(Y,size(F,1),3);

end

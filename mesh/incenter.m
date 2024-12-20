function [I,B,R] = incenter(V,F)
  % INCENTER Compute the incenter of a triangle mesh
  %
  % I = incenter(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   I  #F by 3 list of incenters
  %   B  #F by 3 list of barycentric coordinates of incenter
  %   R  #F by 1 list of inradii
  %
  % See also: barycenter, orthocenter, circumradius
  %

  l = edge_lengths(V,F);
  B = l./sum(l,2);

  if nargout >2
    A = doublearea(V,F);
    % perimeter
    S = sum(l,2);
    R = (A./S);
  end

  I = B(:,1).*V(F(:,1),:) + B(:,2).*V(F(:,2),:) + B(:,3).*V(F(:,3),:);
end

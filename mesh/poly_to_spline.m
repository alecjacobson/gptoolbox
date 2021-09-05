function [P,C] = spline_to_poly(V,E)
  % SPLINE_TO_POLY Trivially convert a poly-line into cubic Bézier curves.
  % 
  % [P,C] = spline_to_poly(V,E)
  % 
  % Inputs:
  %   V  #V by dim list of poly-line vertex positions
  %   E  #E by 2 list of edge indices into rows of V
  % Outputs:
  %   P  #V+2*#E by dim list of cubic bezier curve control points 
  %   C  #E by 4 list of control point indices into rows of P
  P = [V;
    V(E(:,1),:)*(2/3) + V(E(:,2),:)*(1/3); ...
    V(E(:,1),:)*(1/3) + V(E(:,2),:)*(2/3)];
  C = [E(:,1) size(V,1)+reshape((1:2*size(E,1)),[],2) E(:,2)];
end

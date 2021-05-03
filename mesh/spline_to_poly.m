function [V,E,I] = spline_to_poly(P,C,tol)
  % SPLINE_TO_POLY Evaluate a cubic Bezier spline as a polyline where each
  % segment corresponds to a locally flat segment of the curve up to given
  % tolerance
  % 
  % [V,E] = spline_to_poly(P,C,tol)
  % 
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  %   tol  tolerance 
  % Outputs:
  %   V  #V by dim list of vertex locations
  %   E  #E by dim list of edge indices into V
  %   I  #E list of indices into 1:#C
  %

  V = P;
  E = [];
  I = [];
  % consider each cubic
  for c = 1:size(C,1)
    Pc = cubic_flat_eval(P(C(c,:),:),tol);
    J = [C(c,1) size(V,1)+(1:size(Pc,1)-2) C(c,4)];
    Ec = [J(1:end-1);J(2:end)]';
    E = [E; Ec];
    V = [V;Pc(2:end-1,:)];
    I = [I;repmat(c,size(Ec,1),1)];
  end
  [V,~,~,E] = remove_unreferenced(V,E);

end

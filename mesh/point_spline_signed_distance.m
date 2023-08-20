function [D,I,S,K] = point_spline_signed_distance(Q,P,C,tol)
  [sqrD,I,S,K] = point_spline_squared_distance(Q,P,C,tol);
  W = spline_winding_number(P,C,Q);
  D = (1-abs(W)*2).*sqrt(sqrD);
end

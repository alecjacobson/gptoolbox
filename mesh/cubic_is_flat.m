function flag = cubic_is_flat(C,tol)
  % CUBIC_IS_FLAT "Piecewise Linear Approximation of Bézier Curves" [Fischer
  % 2000]
  %
  % flag = cubic_is_flat(C,tol)
  %
  % Inputs:
  %   C  4 by dim list of control points
  %   tol  tolerance
  % Outputs:
  %   flag  whether cubic is effectively flat
  %
  % See also: cubic_eval

  ux = (3.0*C(2,1) - 2.0*C(1,1) - C(4,1)).^2;
  uy = (3.0*C(2,2) - 2.0*C(1,2) - C(4,2)).^2;
  vx = (3.0*C(3,1) - 2.0*C(4,1) - C(1,1)).^2;
  vy = (3.0*C(3,2) - 2.0*C(4,2) - C(1,2)).^2;
  ux(ux < vx) = vx(ux < vx);
  uy(uy < vy) = vy(uy < vy);
  tolerance = 16*tol^2;
  flag = (ux+uy) <= tolerance;
end

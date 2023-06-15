function C = quadratic_to_cubic(Q)
  % C = quadratic_to_cubic(Q)
  %
  % Inputs:
  %   Q  3 by dim list of quadratic Bezier control points
  % Ouputs:
  %   C  4 by dim list of cubicBezier control points
  %

  C(1, :) = Q(1, :);
  C(2, :) = (2/3)*Q(1, :) + (1/3)*Q(2, :);
  C(3, :) = (1/3)*Q(2, :) + (2/3)*Q(3, :);
  C(4, :) = Q(3, :);

end

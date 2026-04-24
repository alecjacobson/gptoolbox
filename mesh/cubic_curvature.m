function [k,KN] = cubic_curvature(C,t)
  % Inputs:
  %   4 by dim matrix of control points
  %   t #t by 1 vector of parameter values
  % Outputs:
  %   k #t by 1 vector of curvature values
  %   KN #t by dim matrix of curvature normal vectors
  %
   dim = size(C,2);
  t = t(:);

  % First derivative control polygon (quadratic)
  D1 = 3*(C(2:4,:) - C(1:3,:));      % 3x dim

  % Second derivative control polygon (linear)
  D2 = 6*(C(3:4,:) - 2*C(2:3,:) + C(1:2,:));  % 2x dim

  % Evaluate derivatives using Bernstein basis
  u = 1 - t;

  % r'(t)
  r1 = ...
      (u.^2)      * D1(1,:) + ...
      (2*u.*t)    * D1(2,:) + ...
      (t.^2)      * D1(3,:);

  % r''(t)
  r2 = ...
      u * D2(1,:) + ...
      t * D2(2,:);

  % Speed
  speed = sqrt(sum(r1.^2,2));

  % Avoid divide-by-zero
  eps_mask = speed > 0;

  k  = zeros(size(t));
  KN = zeros(numel(t),dim);

  if dim == 2
    % 2D curvature (scalar cross product)
    cross_val = r1(:,1).*r2(:,2) - r1(:,2).*r2(:,1);
    k(eps_mask) = cross_val(eps_mask) ./ speed(eps_mask).^3;

  elseif dim == 3
    cross_val = cross(r1,r2,2);
    k(eps_mask) = sqrt(sum(cross_val(eps_mask,:).^2,2)) ./ speed(eps_mask).^3;

  else
    % General dimension: project r2 onto normal plane
    T = zeros(size(r1));
    T(eps_mask,:) = r1(eps_mask,:) ./ speed(eps_mask);

    r2_normal = r2 - sum(r2.*T,2).*T;
    k(eps_mask) = sqrt(sum(r2_normal(eps_mask,:).^2,2)) ./ speed(eps_mask).^2;
  end

  % Compute curvature normal vector
  T = zeros(size(r1));
  T(eps_mask,:) = r1(eps_mask,:) ./ speed(eps_mask);

  r2_normal = r2 - sum(r2.*T,2).*T;
  norm_r2n = sqrt(sum(r2_normal.^2,2));

  mask2 = norm_r2n > 0;
  N = zeros(size(r1));
  if any(mask2)
    N(mask2,:) = r2_normal(mask2,:) ./ norm_r2n(mask2);
  end

  KN = N .* abs(k);

end

function [I,T] = spline_inflection_points(P,C)
  % [I,T] = spline_inflection_points(P,C)
  %
  % Inputs:
  %    P: Nx2 array of control points
  %    C: Mx4 array of control point indices for each cubic Bezier curve
  % Outputs:
  %    I: Kx1 array of curve indices for each inflection point
  %    T: Kx1 array of t values for each inflection point

  function T = quadratic_roots(coefs)
    % Solves a*t^2 + b*t + c = 0
    a = coefs(:,3);
    b = coefs(:,2);
    c = coefs(:,1);
    D = b.^2 - 4*a.*c;
    T = NaN(size(coefs,1),2);
    real_mask = D >= 0 & a ~= 0;
    sqrtD = sqrt(D(real_mask));
    T(real_mask,1) = (-b(real_mask) - sqrtD) ./ (2*a(real_mask));
    T(real_mask,2) = (-b(real_mask) + sqrtD) ./ (2*a(real_mask));
    T = sort(T,2);
  end

  % https://stackoverflow.com/a/35906870/148668


  % Given a cubic Bezier curve with control points C1, C2, C3, C4,
  x1 = P(C(:,1),:);
  x2 = 3*(P(C(:,2),:) - P(C(:,1),:));
  x3 = 3*(P(C(:,3),:) - 2*P(C(:,2),:) + P(C(:,1),:));
  x4 = P(C(:,4),:) - 3*P(C(:,3),:) + 3*P(C(:,2),:) - P(C(:,1),:);
  % The inflection point is when the curvature is zero. So in 2D,
  %
  % X(t)' Y(t)'' - Y(t)' X(t)'' = 0
  % a*t^2 + b*t + c = 0
  % a := 3 * (cf[2].X *cf[3].Y - cf[2].Y *cf[3].X);
  % b := 3 * (cf[1].X *cf[3].Y - cf[1].Y *cf[3].X);
  % c := cf[1].X *cf[2].Y - cf[1].Y *cf[2].X;
  % a = 3 .* (x3(:,1) .* x4(:,2) - x3(:,2) .* x4(:,1));
  % b = 3 .* (x2(:,1) .* x4(:,2) - x2(:,2) .* x4(:,1));
  % c = x2(:,1) .* x3(:,2) - x2(:,2) .* x3(:,1);
  coefs = [ ...
   3 .* (x3(:,1) .* x4(:,2) - x3(:,2) .* x4(:,1)), ...
   3 .* (x2(:,1) .* x4(:,2) - x2(:,2) .* x4(:,1)), ...
   x2(:,1) .* x3(:,2) - x2(:,2) .* x3(:,1)];
  T = fast_roots(coefs,0,1);
  % Sort them by occurance along curves and then by curve index
  T = sort(T,2);
  [J,I] = find(~isnan(T'));
  I = reshape(I,[],1);
  T = reshape(T(sub2ind(size(T),I,J)),[],1);

end

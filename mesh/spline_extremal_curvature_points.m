function [I,T,F] = spline_extremal_curvature_points(P,C)
  % [I,T,F] = spline_extremal_curvature_points(P,C)
  %
  % Inputs:
  %    P: Nx2 array of control points
  %    C: Mx4 array of control point indices for each cubic Bezier curve
  % Outputs:
  %    I: Kx1 array of curve indices for each point
  %    T: Kx1 array of t values for each point
  %    F: Kx1 flag: minimal (-1), maximal (+1), inflection (0)
  %
  % See also: spline_inflection_points
  %

  % Cubic power basis:
  % x(t) = x1 + x2 t + x3 t^2 + x4 t^3
  x1 = P(C(:,1),:);
  x2 = 3*(P(C(:,2),:) - P(C(:,1),:));
  x3 = 3*(P(C(:,3),:) - 2*P(C(:,2),:) + P(C(:,1),:));
  x4 = P(C(:,4),:) - 3*P(C(:,3),:) + 3*P(C(:,2),:) - P(C(:,1),:);

  % A(t) = det(x'(t),x''(t))
  A = [ ...
    3 .* (x3(:,1).*x4(:,2) - x3(:,2).*x4(:,1)), ...
    3 .* (x2(:,1).*x4(:,2) - x2(:,2).*x4(:,1)), ...
         (x2(:,1).*x3(:,2) - x2(:,2).*x3(:,1))];

  % B(t) = ||x'(t)||^2
  B = [ ...
    9  .* sum(x4.*x4,2), ...
    12 .* sum(x3.*x4,2), ...
    4  .* sum(x3.*x3,2) + 6 .* sum(x2.*x4,2), ...
    4  .* sum(x2.*x3,2), ...
          sum(x2.*x2,2)];

  % Inflections: A(t) = 0
  T_inf = fast_roots(A,0,1);

  % Stationary signed curvature:
  %   d/dt A/B^(3/2) = 0
  % gives:
  %   2 A' B - 3 A B' = 0
  Ap = polyder_mat(A);
  Bp = polyder_mat(B);
  Q = 2*conv_mat(Ap,B) - 3*conv_mat(A,Bp);

  T_ext = fast_roots(Q,0,1);

  % Classify curvature extrema
  eps_t = 1e-6;
  F_ext = nan(size(T_ext));

  for i = 1:size(T_ext,1)
    for j = 1:size(T_ext,2)
      t = T_ext(i,j);
      if isnan(t)
        continue;
      end

      tl = max(0,t-eps_t);
      tr = min(1,t+eps_t);

      k0 = signed_curvature_power(x2(i,:),x3(i,:),x4(i,:),t);
      kl = signed_curvature_power(x2(i,:),x3(i,:),x4(i,:),tl);
      kr = signed_curvature_power(x2(i,:),x3(i,:),x4(i,:),tr);

      if kl > k0 && kr > k0
        F_ext(i,j) = -1; % local minimum
      elseif kl < k0 && kr < k0
        F_ext(i,j) = +1; % local maximum
      end
    end
  end

  T_ext(isnan(F_ext)) = nan;

  % Inflection flags
  F_inf = zeros(size(T_inf));

  % Combine extrema and inflections
  T_all = [T_ext T_inf];
  F_all = [F_ext F_inf];

  % Sort within each curve by t
  [T_all,perm] = sort(T_all,2);

  rows = repmat((1:size(T_all,1))',1,size(T_all,2));
  F_all = F_all(sub2ind(size(F_all),rows,perm));

  % Flatten, sorted by curve index then occurrence along curve
  [J,I] = find(~isnan(T_all'));
  I = reshape(I,[],1);

  ind = sub2ind(size(T_all),I,J);
  T = reshape(T_all(ind),[],1);
  F = reshape(F_all(ind),[],1);
end

function D = polyder_mat(C)
  deg = size(C,2)-1;
  D = C(:,1:end-1) .* (deg:-1:1);
end

function C = conv_mat(A,B)
  n = size(A,1);
  p = size(A,2);
  q = size(B,2);

  C = zeros(n,p+q-1);

  for k = 1:p
    C(:,k:k+q-1) = C(:,k:k+q-1) + A(:,k).*B;
  end
end

function k = signed_curvature_power(x2,x3,x4,t)
  v = x2 + 2*x3*t + 3*x4*t.^2;
  a = 2*x3 + 6*x4*t;

  num = v(1).*a(2) - v(2).*a(1);
  den = sum(v.^2).^(3/2);

  k = num ./ den;
end

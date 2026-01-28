function [flag,T] = segment_segment_intersection(A,B,C,D)
  % [flag,T] = segment_segment_intersection(A,B,C,D)
  %
  % Inputs:
  %   A  #A by 2 list of first segment startpoints
  %   B  #A by 2 list of first segment endpoints
  %   C  #C by 2 list of second segment startpoints
  %   D  #C by 2 list of second segment endpoints
  % Outputs:
  %   flag  max(#A,#B) flags of whether segments intersect
  %   T  #flag by 2 list of intersection parameters so that the intersection
  %     occurs at A+T(:,1).*(B-A) and C+T(:,2).*(D-C).

  function flag = check_colinear(o1,A,B,C)
    flag = false(size(o1));
    maybe = o1 == 0;
    % min(ax,bx) ≤ cx ≤ max(ax,bx)  &&  min(ay,by) ≤ cy ≤ max(ay,by)
    flag(maybe) = ...
      min(A(maybe,1),B(maybe,1)) <= C(maybe,1) & ...
      C(maybe,1) <= max(A(maybe,1),B(maybe,1)) & ...
      min(A(maybe,2),B(maybe,2)) <= C(maybe,2) & ...
      C(maybe,2) <= max(A(maybe,2),B(maybe,2));
  end
  function z = cross2(X, Y)
    z = X(:,1).*Y(:,2) - X(:,2).*Y(:,1);
  end
  function T = proper_params(A,B,C,D)
    % Assume we already know that (A,B) intersects (C,D) properly.
    %
    % Inputs:
    %   A  #A by 2 list of first segment startpoints
    %   B  #A by 2 list of first segment endpoints
    %   C  #C by 2 list of second segment startpoints
    %   D  #C by 2 list of second segment endpoints
    % Outputs:
    %   T  #T by 2 list of intersection parameters so that the intersection
    %     occurs at A+T(:,1).*(B-A) and C+T(:,2).*(D-C).
    %
    r = B - A;
    s = D - C;
    qmp = C - A;
    denom = cross2(r, s);
    t = cross2(qmp, s) ./ denom;
    u = cross2(qmp, r) ./ denom;
    T = [t u];
  end
  function T = colinear_param(A,B,C)
    % Assume we already know that C lies on (A,B).
    %
    % Inputs:
    %   A  #A by 2 list of first segment startpoints
    %   B  #A by 2 list of first segment endpoints
    %   C  #C by 2 list of second segment startpoint 
    % Outputs:
    %   T  #T by 1 so that C = A + T.*(B-A).
    T = zeros(size(C,1),1);
    % T = (C-A)⋅(B-A) / ||B-A||²
    r = B - A;
    T = sum((C - A) .* r, 2) ./ sum(r.^2, 2);
  end

  o1 = orient2d(A,B,C);
  o2 = orient2d(A,B,D);
  o3 = orient2d(C,D,A);
  o4 = orient2d(C,D,B);
  proper = (o1 > 0 & o2 < 0 | o1 < 0 & o2 > 0) & (o3 > 0 & o4 < 0 | o3 < 0 & o4 > 0);

  colinear_1 = check_colinear(o1,A,B,C);
  colinear_2 = check_colinear(o2,A,B,D);
  colinear_3 = check_colinear(o3,C,D,A);
  colinear_4 = check_colinear(o4,C,D,B);

  flag = proper | colinear_1 | colinear_2 | colinear_3 | colinear_4;
  if nargout
    T = nan(numel(flag),2);
    T(proper,:) = proper_params(A(proper,:),B(proper,:),C(proper,:),D(proper,:));
    T(colinear_1,1) = colinear_param(A(colinear_1,:),B(colinear_1,:),C(colinear_1,:));
    T(colinear_1,2) = 0;

    T(colinear_2,1) = colinear_param(A(colinear_2,:),B(colinear_2,:),D(colinear_2,:));
    T(colinear_2,2) = 1;

    T(colinear_3,2) = colinear_param(C(colinear_3,:),D(colinear_3,:),A(colinear_3,:));
    T(colinear_3,1) = 0;

    T(colinear_4,2) = colinear_param(C(colinear_4,:),D(colinear_4,:),B(colinear_4,:));
    T(colinear_4,1) = 1;

  end

end

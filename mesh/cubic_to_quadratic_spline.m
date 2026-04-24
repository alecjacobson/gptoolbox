function [P,Q] = cubic_to_quadratic_spline(C)
  % [P,Q] = cubic_to_quadratic_spline(C)

  method = 'naive';
  switch method
  case 'naive'
  case 'minimal' % not a good name
    [I,T] = spline_inflection_points(C,[1 2 3 4]);
    switch numel(T)
    case {0,1}
      switch numel(T)
      case 0
        % t=0.5 seems hard to beat significantly
        t = 0.5;
        f1 = 0.5;
        f2 = 0.5;
        M = cubic_eval(C,t);
        V = cubic_tangent(C,t);
        [~,valid1,f1] = intersect_lines_2d(C(1,:),C(2,:)-C(1,:),M,V);
        [~,valid2,f2] = intersect_lines_2d(C(4,:),C(3,:)-C(4,:),M,V);
      case 1
        t = T;
        % Not sure why, but this always seems to work at inflection points. No need
        % to intersect_lines_2d.
        f1 = t;
        f2 = 1-t;
        M = cubic_eval(C,t);
      end
      Q1 = [C(1,:);(C(2,:)-C(1,:))*f1 + C(1,:);M];
      Q2 = [M;     (C(3,:)-C(4,:))*f2 + C(4,:);C(4,:)];
      P = [Q1;Q2(2:end,:)];
      Q = [1 2 3;3 4 5];
    case 2
      t1 = T(1);
      t2 = T(2);
      M1 = cubic_eval(C,t1);
      V1 = cubic_tangent(C,t1);
      M2 = cubic_eval(C,t2);
      V2 = cubic_tangent(C,t2);

      [~,valid1,f1] = intersect_lines_2d(C(1,:),C(2,:)-C(1,:),M1,V1);
      [~,valid2,f2] = intersect_lines_2d(M1,V1,M2,V2);
      [~,valid3,f3] = intersect_lines_2d(C(4,:),C(3,:)-C(4,:),M2,V2);
      Q1 = [C(1,:);(C(2,:)-C(1,:))*f1 + C(1,:);M1];
      Q2 = [M1;V1*f2 + M1;M2];
      Q3 = [M2;     (C(3,:)-C(4,:))*f3 + C(4,:);C(4,:)];
      P = [Q1;Q2(2:end,:);Q3(2:end,:)];
      Q = [1 2 3;3 4 5;5 6 7];
    end
  end
end

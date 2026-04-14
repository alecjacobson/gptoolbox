function [T0,T1] = spline_robust_endpoint_tangents(P,C)
  % [T0,T1] = spline_robust_endpoint_tangents(P,C)
  function T0 = helper_1(C1,C4)
    T0 = C4 - C1;
  end
  function T0 = helper_2(C1,C3,C4,epsilon)
    T0 = C3 - C1;
    too_small = sum(T0.^2,2) < epsilon;
    T0(too_small,:) = helper_1(C1(too_small,:),C4(too_small,:));
  end
  
  function T0 = helper_3(C1,C2,C3,C4,epsilon)
    T0 = C2 - C1;
    too_small = sum(T0.^2,2) < epsilon;
    T0(too_small,:) = helper_2(C1(too_small,:),C3(too_small,:),C4(too_small,:),epsilon);
  end
  epsilon = 1e-7;
  C1 = P(C(:,1),:);
  C2 = P(C(:,2),:);
  C3 = P(C(:,3),:);
  C4 = P(C(:,4),:);
  T0 = helper_3(C1,C2,C3,C4,epsilon);
  T1 = -helper_3(C4,C3,C2,C1,epsilon);
end

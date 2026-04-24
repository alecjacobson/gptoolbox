function [k14] = spline_endpoint_curvature(P,C)
  A123 = doublearea(P,C(:,[1 2 3]));
  A234 = doublearea(P,C(:,[2 3 4]));
  s1 = normrow(P(C(:,2),:) - P(C(:,1),:)).^3;
  s4 = normrow(P(C(:,4),:) - P(C(:,3),:)).^3;
  k14 = (4/6).*[A123./s1, A234./s4];
end

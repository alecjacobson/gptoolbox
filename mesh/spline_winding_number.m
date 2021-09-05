function W = spline_winding_number(P,C,O)
  % SPLINE_WINDING_NUMBER Compute generalized winding numbers for a query
  % points O with respect to a given cubic Bezier spline (P,C).
  %
  % W = spline_winding_number(P,C,O)
  %
  % Inputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  %   O  #O by 2 list of query points
  % Outputs:
  %   W  #O list of winding number values
  %

  %% Unnecessary? This is not a bottleneck of the flattening
  %% Is entire spline pwn ?
  %E = C(:,[1 4]);
  %pwn = all(full(sparse(E,1,repmat([1 -1],size(E,1),1))));
  %if pwn
  %  % then immediately identify everything outside of the full convex hull
  %  H = convhull(Ci);
  %  I = find(~inpolygon(V(:,1),V(:,2),C(H(1:end-1),1),C(H(1:end-1),2)));
  %else
  %  I = 1:size(V,1);
  %end
  %W = zeros(size(O,1),1);
  %I = 1:size(V,1);
  %for ci = 1:size(C,1)
  %  Ci = P(C(ci,:),:);
  %  W(I) = W(I) + cubic_winding_number(Ci,O(I,:));
  %end

  W = zeros(size(O,1),1);
  for ci = 1:size(C,1)
    Ci = P(C(ci,:),:);
    W = W + cubic_winding_number(Ci,O);
  end

end

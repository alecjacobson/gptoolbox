function [sqrD,I,S,K,T] = point_spline_squared_distance(Q,P,C,tol)
  % [sqrD,I,S,K,T] = point_spline_squared_distance(Q,P,C,tol)
  %
  % Inputs:
  %   Q   #Q by dim list of query points
  %   P   #P by dim list of control points
  %   C   #C by 4 list of indices into P of cubic splines
  %   tol  tolerance for is_flat
  % Outputs:
  %   sqrD  #Q list of squared distances
  %   I  #Q list of indices into C
  %   S  #Q list of parameter value on corresponding cubic spline
  %   K  #Q by dim list of closest points on corresponding cubic spline
  %   T  #Q by dim list of tangent vectors on corresponding cubic spline
  %   
  % O(nm)
  sqrD = inf(size(Q,1),1);
  I = nan(size(Q,1),1);
  S = nan(size(Q,1),1);
  K = nan(size(Q));
  T = nan(size(Q));
  H = [];

  % This is still O(nm)...
  % max(d_i) ≥ d
  for c = 1:size(C,1)
    Bc = pdist2(Q,P(C(c,:),:),'squaredeuclidean');
    sqrD = max(sqrD,max(Bc,[],2));
  end
  high = sqrD;


  % Don't even compute the convex hull, just immediately compute distance to all
  % 6 possible edges
  E = nchoosek(1:4,2);
  for c = 1:size(C,1)
    % d_hull ≤ d
    Pc = P(C(c,:),:);
    low = point_mesh_squared_distance(Q,Pc,E);
    J = find(low < sqrD);

    if ~isempty(J)
      [sqrDc,Sc] = point_cubic_squared_distance(Q(J,:),Pc,tol);
      smaller = sqrDc < sqrD(J);
      Jsmaller = J(smaller);
      sqrD(Jsmaller) = sqrDc(smaller);
         I(Jsmaller) = c;
         S(Jsmaller) = Sc(smaller);
      if nargout>3
        K(Jsmaller,:) = cubic_eval(Pc,S(Jsmaller));
      end
      if nargout>4
        T(Jsmaller,:) = cubic_tangent(Pc,S(Jsmaller));
      end
    end
  end
end

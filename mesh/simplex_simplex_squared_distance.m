function [sqrd,ca,cb] = simplex_simplex_squared_distance(A,B)
  % SIMPLEX_SIMPLEX_SQUARED_DISTANCE Squared distance between two simplices.
  %
  % [sqrd,ca,cb] = simplex_simplex_squared_distance(A,B)
  %
  % Inputs:
  %   A  #A by dim list of vertices of a (#A-1)-simplex
  %   B  #B by dim list of vertices of a (#B-1)-simplex
  % Outputs:
  %   sqrd  squared distance between the simplices
  %   ca  1 by dim closest point on A
  %   cb  1 by dim closest point on B
  %
  % The simplices may have different dimensions. Degenerate simplices are
  % also handled.
  %

  assert(size(A,2) == size(B,2));
  assert(~isempty(A));
  assert(~isempty(B));

  sqrd = inf;
  ca = [];
  cb = [];

  recursive_helper(A,B);

  function recursive_helper(A,B)

    % Nothing can improve on an intersection.
    if sqrd == 0
      return;
    end

    na = size(A,1);
    nb = size(B,1);

    % Point-point.
    if na == 1 && nb == 1
      d = sum((A-B).^2);
      if d < sqrd
        sqrd = d;
        ca = A;
        cb = B;
      end
      return;
    end

    % Parameterize the two affine hulls:
    %
    %   a = A(1,:) + u'*EA
    %   b = B(1,:) + v'*EB
    %
    % and find their closest pair by least squares.
    EA = A(2:end,:) - A(1,:);
    EB = B(2:end,:) - B(1,:);

    M = [EA' -EB'];
    x = pinv(M) * (B(1,:)-A(1,:))';

    u = x(1:na-1);
    v = x(na:end);

    wa = [1-sum(u);u];
    wb = [1-sum(v);v];

    candidate_ca = wa'*A;
    candidate_cb = wb'*B;
    affine_sqrd = sum((candidate_ca-candidate_cb).^2);

    % Distance between the affine hulls is a lower bound for this entire
    % recursive problem.
    if affine_sqrd >= sqrd
      return;
    end

    % If the affine-hull closest pair lies in both simplices, then it is also
    % the closest pair of the simplices.
    tol = 1e-12;
    if all(wa >= -tol) && all(wb >= -tol)
      sqrd = affine_sqrd;
      ca = candidate_ca;
      cb = candidate_cb;
      return;
    end

    % Otherwise a closest pair occurs on the boundary of at least one
    % simplex. Recurse over all codimension-one facets.
    if na > 1
      for i = 1:na
        I = [1:i-1 i+1:na];
        recursive_helper(A(I,:),B);
        if sqrd == 0
          return;
        end
      end
    end

    if nb > 1
      for i = 1:nb
        I = [1:i-1 i+1:nb];
        recursive_helper(A,B(I,:));
        if sqrd == 0
          return;
        end
      end
    end
  end
end

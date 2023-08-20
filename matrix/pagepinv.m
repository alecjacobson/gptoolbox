function Y = pagepinv(X,tol)
  % Y = pagepinv(X,tol)
  %
  % Inputs:
  %    X  m x n x p array
  % Outputs:
  %    Y  n x m x p array so that Y(:,:,i) = pinv(X(:,:,i))
  %
  [U,s,V] = pagesvd(X,'econ','vector');

  if nargin<2
    tol = max(size(X(:,:,1))) * eps(pagenorm(s,inf));
  end

  s(s>tol) = 1./s(s>tol);

  Y = pagemtimes(permute(s,[2 1 3]).*V,'none',U,'ctranspose');
end

function Y = pagenull(X,tol)
  [~,s,V] = pagesvd(X,'vector');

  if nargin < 2
    tol = max(size(X,1),size(X,2)) .* eps(max(s,[],1));
  end

  r = sum(s > tol,1);

  if any(r ~= r(1))
    error('pagenull:variableNullity', ...
      'Pages have different null-space dimensions.');
  end

  Y = V(:,r(1)+1:end,:);
end
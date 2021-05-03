function H = psd_project_rows(H,tol)
  % H = psd_project_rows(H)
  if nargin < 2
    tol = 0;
  end
  m = round(sqrt(size(H,2)));
  assert(m*m==size(H,2));
  for i = 1:size(H,1)
    A = reshape(H(i,:),m,m);
    [V,D] = eig(triu(A)+triu(A,1)','vector');
    B = V*(max(D,tol).*V');
    H(i,:) = B(:);
    if any(~isreal(B(:)))
      keyboard
    end
  end

end

function S = stiffness_matrix(n)
  i = [1:n, 1:n-1, 2:n];
  j = [1:n, 2:n, 1:n-1];
  v = [[-2,-4*ones(1,n-2),-2], ones(1,n-1), ones(1,n-1)];
  S = sparse(i,j,v,n,n);
end 
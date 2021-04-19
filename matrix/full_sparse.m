function F = full_sparse(I,J,V,S1,S2)
  % FULL_SPARSE Faster than calling:
  %
  % F = full_sparse(I,J,V,S1,S2)
  %
  % F = full(sparse(I,J,V,S1,S2))
  %
  % See also: sparse, accumarray
  if nargin>3
    F = accumarray([I(:) J(:)],V(:),[S1 S2]);
  else
    F = accumarray([I(:) J(:)],V(:));
  end
end

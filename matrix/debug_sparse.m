function S = debug_sparse(I,J,V,varargin)
  fprintf('I: %d %d\n',size(I));
  fprintf('J: %d %d\n',size(J));
  fprintf('V: %d %d\n',size(V));
  if nargin>3
    m = varargin{1};
    n = varargin{2};
    fprintf('S: %d %d\n',m,n);
  else
    fprintf('S: %d %d\n',max(I(:)),max(J(:)));
  end
  S = sparse(I,J,V,varargin{:});
end

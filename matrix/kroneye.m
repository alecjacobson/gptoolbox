function K = kroneye(A,n)
  % KRONEYE Take kronecker product with sparse identity in fast way.  Compare
  % performance to kron(A,speye(n,n))
  %
  % K = kroneye(A,n)
  %
  % Inputs:
  %   A  #rows by #cols sparse matrix
  %   n  size of identity
  % Outputs:
  %   K  #rows*n by #cols*n sparse matrix
  %
  % See also: repdiag, kron
  %

  %% Slow:
  %% http://www.mathworks.com/matlabcentral/newsreader/view_thread/17046
  %K(1:m+1:m^2,:) = repmat(A(:).',n,1);
  %K = reshape(permute(reshape(K,[m n size(A)]),[1 3 2 4]),[m*size(A,1) n*size(A,2)]);
  if ~issparse(A)
    A = sparse(A);
  end
  
  if n >10*size(A,1)
    K = kron(A,speye(n));
    return;
  end

  % Compute kron(speye(n,n),A) in a fast way
  C = cell(n,1);
  [C{:}] = deal(A);
  K = blkdiag(C{:});
  % Rearrange rows and columns
  I = reshape(reshape(1:(size(A,1)*n),size(A,1),n)',size(A,1)*n,1);
  J = reshape(reshape(1:(size(A,2)*n),size(A,2),n)',size(A,2)*n,1);
  K = K(I,J);
end

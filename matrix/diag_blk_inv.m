function N = diag_blk_inv(M,k)
  % DIAG_BLK_INV Invert a block matrix made out of square diagonal blocks. M
  % should be of the form:
  %
  %   M = [M11 M12 ... M1n
  %        M21 ...
  %        ...
  %        Mn1 Mn2 ... Mnn];
  % where each Mij is a diagonal matrix.
  %
  % Inputs:
  %   M  n*k by n*k matrix made out of diagonal blocks.
  %   k  size of blocks
  % Outputs:
  %   N  n*k by n*k matrix inverse of M
  %

  switch k
  case 1
    assert(isdiag(M));
    N = diag(sparse(1./diag(M)));
    return;
  end

  assert(size(M,1)==size(M,2),'Must be square');
  assert(rem(size(M,1),k)==0,'Must be divisble by k');
  n = size(M,1)/k;
  % Extract 4 blocks
  A = M(0*n+(1:(k-1)*n),0*n+(1:(k-1)*n));
  B = M(0*n+(1:(k-1)*n),(k-1)*n+(1:n));
  C = M((k-1)*n+(1:n),0*n+(1:(k-1)*n));
  D = M((k-1)*n+(1:n),(k-1)*n+(1:n));
  assert(isdiag(D));
  % https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion
  Ainv = diag_blk_inv(A,k-1);
  % Schur complement 
  S = (D-C*Ainv*B);
  assert(isdiag(S));
  Sinv = diag_blk_inv(S,1);
  N = [Ainv + Ainv*B*Sinv*C*Ainv -Ainv*B*Sinv; ...
       -Sinv*C*Ainv              Sinv];
end

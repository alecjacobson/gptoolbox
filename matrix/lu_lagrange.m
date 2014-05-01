function [L,U,p,PT] = lu_lagrange(ATA,C,J,S)
  % LU_LAGRANGE Compute a LU decomposition for a special type of
  % matrix Q that is symmetric but not positive-definite:
  % Q = [A'*A C
  %      C'   0];
  % where A'*A, or ATA, is given as a symmetric positive definite matrix and C
  % has full column-rank(?)
  % http://www.alecjacobson.com/weblog/?p=2242
  % http://www.alecjacobson.com/weblog/?p=2480
  %
  % This is often called the "Schur complement trick".
  %
  % [L,U] = lu_lagrange(ATA,C)
  % [L,U,p,PT] = lu_lagrange(ATA,C,J,S)
  %
  % Inputs:
  %   ATA   n by n square, symmetric, positive-definite system matrix, usually
  %     the quadratic coefficients corresponding to the original unknowns in a
  %     system
  %   C  n by m rectangular matrix corresponding the quadratic coefficients of
  %     the original unknowns times the lagrange multipliers enforcing linear
  %     equality constraints
  %   Optional:
  %     J  n by n lower triangular result of [J] = chol(ATA,'lower')
  %     S  n by n permutation matrix result of [J,p,S] = chol(ATA,'lower'),
  %       notice now J is also different, requires that PT is computed in
  %       output parameters
  % Outputs:
  %   L  lower triangular matrix such that Q = L*U
  %   U  upper triangular matrix such that Q = L*U
  %   p  flag returned from last call to chol
  %   PT  permutation matrix, such that PT'*Q*PT = L*U, requires that S is given
  %     if J is given
  %
  % Notes: This only seems worthwhile when m << n
  %
  % See also: chol, lu, min_quad_with_fixed
  % 

  % number of unknowns
  n = size(ATA,1);
  % number of lagrange multipliers
  m = size(C,2);
  % constraints should not have null columns
  assert(~any(all(C==0)));

  assert(size(ATA,2) == n);
  assert(size(C,1) == n);

  L = [];
  U = [];
  p = 0;
  PT = [];

  if ~exist('J','var') || isempty(J)
    if nargout >= 4
      % compute lower triangular of cholesky factorization of ATA, with
      % permutation matrix
      [J,p,S] = chol(ATA,'lower');
    else
      % compute lower triangular of cholesky factorization of ATA
      [J,p] = chol(ATA,'lower');
    end
  else
    assert(all(size(J) == size(ATA)))
    % if permutation matrix S of chol factor of A is given then we must be
    % outputing final permutation S
    if exist('S','var')
      assert(nargout >= 4)
    end
  end

  % if computing permutation matrix in output then we need S permutation of
  % chol factor of ATA
  if nargout >= 4
    assert(0 ~= exist('S','var'));
    assert(all(size(S) == size(ATA)));
  end

  if p ~= 0
    return;
  end

  % construct helper matrix M
  if nargout >= 4
    M = J\(S'*C);
  else
    M = J\C;
  end

  % compute cholesky factorization of M'*M
  if nargout >= 4
    if issparse(M)
      [K,p,T] = chol(M'*M,'lower');
    else
      [K,p] = chol(M'*M,'lower');
      T = speye(size(K));
    end
  else
    [K,p] = chol(M'*M,'lower');
  end

  if p ~= 0
    return;
  end

  % assemble LU decomposition of Q
  if(issparse(ATA))
    Z = sparse(n,m);
  else
    Z = zeros(n,m);
  end

  if nargout >= 4
    M = M*T;
    PT = blkdiag(S,T);
  end

  U = [J' M;Z' -K'];
  %L = [J Z;M' K];
  % this is slightly faster
  L = U';
  L(n+(1:m),n+(1:m)) = -L(n+(1:m),n+(1:m));

end

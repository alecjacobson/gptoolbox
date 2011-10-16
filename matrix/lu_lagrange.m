function [L,U,p] = lu_lagrange(ATA,C,J)
  % LU_LAGRANGE Compute a LU decomposition for a special type of
  % matrix Q that is symmetric but not positive-definite:
  % Q = [A'*A C
  %      C'   0];
  % where A'*A, or ATA, is given as a symmetric positive definite matrix and C
  % has full column-rank(?)
  %
  % [J] = lu_lagrange(ATA,C)
  %
  % Inputs:
  %   ATA   n by n square, symmetric, positive-definite system matrix, usually
  %     the quadratic coefficients corresponding to the original unknowns in a
  %     system
  %   C  n by m rectangular matrix corresponding the quadratic coefficients of
  %     the original unknowns times the lagrange multipliers enforcing linear
  %     equality constraints
  %   Optional:
  %     J  n by n lower triangular result of [J] = chol(ATA)
  % Outputs:
  %   L  lower triangular matrix such that Q = L*U
  %   U  upper triangular matrix such that Q = L*U
  %   p  flag returned from last call to chol
  %
  % See also: chol, lu, min_quad_with_fixed
  % 

  L = [];
  U = [];
  p = 0;

  if ~exist('J','var') || isempty(J)
    % compute lower triangular of cholesky factorization of ATA
    [J,p] = chol(ATA,'lower');
  else
    assert(all(size(J) == size(ATA)))
  end

  if p ~= 0
    return;
  end

  % construct helper matrix M
  M = J\C;

  % compute cholesky factorization of M'*M
  [K,p] = chol(M'*M,'lower');

  if p ~= 0
    return;
  end

  % assembly lower triangular of cholesky factorization of Q
  Z = zeros(size(ATA,1),size(C,2));
  L = [J Z;M' K];
  U = [J' M;Z' -K'];
end

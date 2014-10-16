function N = sparseNull(S, tol)
  % spnull returns computes the sparse Null basis of a matrix
  %
  % N = sparseNull(S, tol)
  %
  % Computes a basis of the null space for a sparse matrix.  For sparse
  % matrixes this is much faster than using null.  It does however have lower
  % numerical accuracy.  N is itself sparse and not orthonormal.  So in this
  % way it is like using N = null(S, 'r'), except of course much faster.
  %
  % Jan Schellenberger 10/20/2009
  % based on this:
  % http://www.mathworks.com/matlabcentral/fileexchange/11120-null-space-of-a-sparse-matrix
  if nargin <2
      tol = 1e-9;
  end
  [SpLeft, SpRight] = spspaces(S,2, tol);
  N = SpRight{1}(:,SpRight{3});
  N(abs(N) < tol) = 0;
end

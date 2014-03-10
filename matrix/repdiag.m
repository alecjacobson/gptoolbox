function B = repdiag(A,d)
  % REPDIAG repeat a matrix along the diagonal a certain number of times, so
  % that if A is a m by n matrix and we want to repeat along the diagonal d
  % times, we get a m*d by n*d matrix B such that:
  % B( (k*m+1):(k*m+1+m-1), (k*n+1):(k*n+1+n-1)) = A 
  % for k from 0 to d-1
  %
  % B = repdiag(A,d)
  %
  % Inputs:
  %   A  m by n matrix we are repeating along the diagonal. May be dense or
  %     sparse
  %   d  number of times to repeat A along the diagonal
  % Outputs:
  %   B  m*d by n*d matrix with A repeated d times along the diagonal,
  %     will be dense or sparse to match A
  %
  % See also: kroneye
  %

  %m = size(A,1);
  %n = size(A,2);
  %if(issparse(A))
  %  [I,J,V] = find(A);
  %  BI = I;
  %  BJ = J;
  %  BV = V;
  %  for k = 2:d
  %    BI = [BI (k-1)*m+I];
  %    BJ = [BJ (k-1)*n+J];
  %    BV = [BV V];
  %  end
  %  B = sparse(BI,BJ,BV,m*d,n*d);
  %else
  %  B = zeros(m*d,n*d);
  %  for k = 0:(d-1)
  %    B( (k*m+1):(k*m+1+m-1), (k*n+1):(k*n+1+n-1)) = A;
  %  end
  %end

  % http://www.physicsforums.com/showthread.php?t=77645
  % Also slow:
  % B = kron(speye(d),A);
  % 10x faster than for loop IJV
  C = cell(d,1);
  [C{:}] = deal(A);
  B = blkdiag(C{:});


end

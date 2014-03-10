function SV = stacktimes(S,V)
  % STACKTIMES Multiply each matrix in a vertical stack of matrices, S, times
  % each vector in a list of row-vectors, V. Equivalent to the following:
  %     s = size(V,1);
  %     n = size(V,2);
  %     m = size(S,1)/s;
  %     SV = zeros(s,m);
  %     for ii = 0:(s-1)
  %      SV(ii+1,:) = (S((m*ii+1):(m*ii+m),:)*(V(ii+1,:)'))';
  %     end
  %
  % SV = stacktimes(S,V)
  %
  % Inputs:
  %   S  a list of s, m by n matrices, stacked vertically. That is a s*m by n
  %     matrix
  %      OR
  %      a 3d-array, of size m by n by s
  %   V  a list of s n-size row-vectors, that is: an s by n array
  % Outputs:
  %   SV = list of s m-size row-vectors
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  s = size(V,1);
  n = size(V,2);
  if(size(S,3) ~= 1)
    assert(s == size(S,3));
    m = size(S,1);
    % S is a 3D array, convert to vertically stacked
    S = reshape(permute(S,[2,1,3]),[n,s*m])';
  else
    assert( n == size(S,2));
    m = size(S,1)/s;
    assert( m == round(m));
  end

  I = repmat([1:(m*s)],1,n);
  J = repmat(reshape(repmat(1:n:(s*n),m,1),1,s*m),1,n) + ...
    reshape(repmat(0:(n-1),s*m,1),1,s*m*n);
  SV = reshape(sparse(I,J,S(:))*reshape(V',s*n,1),m,s)';
end

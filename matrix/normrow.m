function [ B ] = normrow( A )
  % NORMROW  Compute l2 row vector norms
  %
  % B = normrow( A )
  %
  % Input:
  %  A  #A by D list of row vectors of dimension D
  % Output:
  %  B  #A list of norms of row vectors in A
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Daniele Panozzo
  %
  B = sqrt(sum(A.^2,2));
end


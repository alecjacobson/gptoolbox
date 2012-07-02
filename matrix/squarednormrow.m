function [ B ] = squarednormrow( A )
  % SQUAREDNORMROW 
  %
  % B = squarednormrow( A )
  %
  % Input:
  %  A  #A by D list of row vectors of dimension D
  % Output:
  %  B  #B list of norms of row vectors in A
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Daniele Panozzo
  %
  B = sum(A.^2,2);
end


function [ A ] = normalizerow( A ) %#codegen
  % NORMALIZEROW Normalize each row so that each row's l2 norm as a vector is 1
  %
  % [ A ] = normalizerow( A )
  %
  % Input:
  %  A  #A by D list of row vectors of dimension D
  % Output:
  %  B  #B by D list of normalized row vectors 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Daniele Panozzo
  %

  if issparse(A)
    % speed up (20x) for large sparse matrices
    A = bsxfun(@times,A,1./sqrt(sum(A.^2,2)));
  else
    A = A./repmat(sqrt(sum(A.^2,2)),1,size(A,2));
  end
end


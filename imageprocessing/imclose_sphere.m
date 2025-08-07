function B = imclose_sphere(A,r)
  % Compute morphological closing with sphere of radius r.
  %
  % B = imclose_sphere(A,r)
  %
  % Inputs:
  %   A  binary image
  %   r  radius of sphere
  % Outputs:
  %   B  binary image same size as A
  % 
  % % Example:
  % M = (imclose(A,strel('sphere',32)));
  % F = imclose_sphere(A,32);
  %
  assert(islogical(A));
  % pad A by r
  rc = ceil(r);
  A = padarray(A,repmat(rc,1,ndims(A)),'both');
  B = bwdist(bwdist(A>0)>r)>r;
  sizeB = size(B);
  indices = arrayfun(@(dim) (rc+1):(sizeB(dim)-rc), 1:ndims(B), 'UniformOutput', false);
  B = B(indices{:});
end

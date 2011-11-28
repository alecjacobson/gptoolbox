function L = wachspress_laplacian(V,F)
  % WACHSPRESS_LAPLACIAN  Discrete laplacian using wachspress weights.
  %
  % L = wachspress_laplacian(V,F)
  %
  % Inputs:
  %   V  #V x dim matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  % Outputs:
  %   L  sparse #V x #V matrix of wachspress weights 
  %
  % See also: cotmatrix, cotmatrix_intrinsic
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Denis Zorin
  %

  % number of vertices
  n = size(V,1);

  % edge lengths numbered same as opposite vertices
  l = [ ...
    sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
    sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
    sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
    ];

  % renaming indices of vertices of triangles for convenience
  i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); 
  l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
  % semiperimeters
  s = (l1 + l2 + l3)*0.5;
  % Heron's formula for area
  dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
  % cotangents and diagonal entries for element matrices
  % Note: same as cotmatrix but divide by edge lengths
  % correctly divided by 4 (alec 2010)
  cot12 = ((l1.^2 + l2.^2 -l3.^2)./dblA/4)./(l3.^2); 
  cot23 = ((l2.^2 + l3.^2 -l1.^2)./dblA/4)./(l1.^2); 
  cot31 = ((l1.^2 + l3.^2 -l2.^2)./dblA/4)./(l2.^2); 
  % diag entries computed from the condition that rows of the matrix sum up to 1
  % (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
  diag1 = (-cot12-cot31); 
  diag2 = (-cot12-cot23); 
  diag3 = (-cot31-cot23);
  % indices of nonzero elements in the matrix for sparse() constructor
  i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
  j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
  % values corresponding to pairs form (i,j)
  v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
  % for repeated indices (i,j) sparse automatically sums up elements, as we want
  L = sparse(i,j,v,n,n);

end


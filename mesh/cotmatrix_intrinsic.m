function L = cotmatrix_intrinsic(l,F,nvert)
  % COTMATRIX_INTRINSIC Compute cotangent Laplacian using edge lengths (not
  % vertex positions).
  %
  % L = cotmatrix(l,F,nvert)
  %
  % Inputs:
  %   l  #F by 3, array of edge lengths of edges opposite each face in F
  %   F  #F by 3, list of indices of triangle corners
  %   nvert  number of vertices, only needed to set size
  % Outputs:
  %   L  sparse nvert x nvert matrix of cot weights 
  %
  % Example:
  %   % Reproduce regular cotmatrix using cotmatrix_intrinsic
  %   % edges numbered same as opposite vertices
  %   FT = F';
  %   l = [ ...
  %     sqrt(sum((V(FT(:,2),:)-V(FT(:,3),:)).^2,2)) ...
  %     sqrt(sum((V(FT(:,3),:)-V(FT(:,1),:)).^2,2)) ...
  %
  % See also: cotmatrix

  % should change code below, so we don't need this transpose
  if(size(F,1) == 3)
    warning('F seems to be 3 by #F, it should be #F by 3');
  end
  F = F';

  % Law of cosines + Law of sine gives you:

  % renaming indices of vertices of triangles for convenience
  i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
  l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
  % semiperimeters
  s = (l1 + l2 + l3)*0.5;
  % Heron's formula for area
  dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
  % cotangents and diagonal entries for element matrices
  % correctly divided by 4 (alec 2010)
  cot12 = (l1.^2 + l2.^2 -l3.^2)./dblA/4; 
  cot23 = (l2.^2 + l3.^2 -l1.^2)./dblA/4; 
  cot31 = (l1.^2 + l3.^2 -l2.^2)./dblA/4; 
  % diag entries computed from the condition that rows of the matrix sum up to 1
  % (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
  diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;
  % indices of nonzero elements in the matrix for sparse() constructor
  i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
  j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
  % values corresponding to pairs form (i,j)
  v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
  % for repeated indices (i,j) sparse automatically sums up elements, as we want
  L = sparse(i,j,v,nvert,nvert);
end

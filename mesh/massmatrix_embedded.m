function M = massmatrix_embedded(V,F, type)
  % MASSMATRIX_EMBEDDED mass matrix for the mesh given by V and F embedded in
  % size(V,2). This is a generalization of massmatrix and a wrapper for
  % massmatrix_intrinsic
  %
  % M = massmatrix_embedded(V,F, type)
  %
  % Inputs:
  %  V  #V x dim matrix of vertex coordinates
  %  F  #F x 3  matrix of indices of triangle corners
  %  type  string containing type of mass matrix to compute
  %   'full': full mass matrix for p.w. linear fem
  %   'barycentric': diagonal lumped mass matrix obtained by summing 1/3
  %   'voronoi': true voronoi area, except in cases where triangle is obtuse
  %     then uses 1/2, 1/4, 1/4
  % Output:
  %  M  #V by #V sparse mass matrix
  %
  % Example:
  %   % (V,F) is a surface mesh in R^3
  %   M = massmatrix(V,F,'voronoi');
  %   Me = massmatrix_embedded(V,F,'voronoi');
  %   % results should be the same up to floating point error
  %   max(max(abs(M-Me)))
  %
  % See also: massmatrix, massmatrix_intrinsic
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

   % edge lengths numbered same as opposite vertices
   l = [ ...
     sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
     sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
     sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
     ];
   M = massmatrix_intrinsic(l,F,size(V,1),type);

end

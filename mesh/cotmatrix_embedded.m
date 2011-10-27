function L = cotmatrix_embedded(V,F)
  % COTMATRIX_EMBEDDED computes cotangent matrix (laplacian mesh operator),
  % (mass/area terms already cancelled out) where (V,F) are a mesh embedded in
  % size(V,2) dimensional space. This is a generalization of cotmatrix and a
  % wrapper for cotmatrix_intrinsic
  %
  % L = cotmatrix_embedded(V,F)
  %
  % Inputs:
  %   V  #V x dim matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  % Outputs:
  %   L  sparse #V x #V matrix of cot weights 
  %
  % Example:
  %   % (V,F) is a surface mesh in R^3
  %   L = cotmatrix(V,F);
  %   Le = cotmatrix_embedded(V,F);
  %   % results should be the same up to floating point error
  %   max(max(abs(L-Le)))
  %
  % See also: cotmatrix, cotmatrix_intrinsic
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Denis Zorin
  %

   % edge lengths numbered same as opposite vertices
   l = [ ...
     sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
     sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
     sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
     ];
   L = cotmatrix_intrinsic(l,F,size(V,1));
end


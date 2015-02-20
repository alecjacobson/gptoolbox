function G = grad3(V,T)
  % GRAD3 Compute the numerical gradient operator for tet meshes in 3d
  % 
  % G = grad(V,F)
  %
  % Inputs:
  %   V  #vertices by dim list of mesh vertex positions
  %   T  #elements by 3 list of mesh tet indices
  % Outputs:
  %   G  #elements*dim by #V Gradient operator
  %
  % Example:
  %   L = cotmatrix(V,T);
  %   G  = grad3(V,T);
  %   vol = volume(V,T);
  %   GMG = -G'*repdiag(diag(sparse(vol)),3)*G;
  %

  warning('Deprecated. Call grad directly.');
  G = grad(V,T);
end

function M = massmatrix3(V,T, type)
  % MASSMATRIX3 mass matrix for the mesh given by V and F
  %
  % M = massmatrix3(V,T, type)
  %
  %
  % Inputs:
  %   V #V x 3 matrix of vertex coordinates
  %   T #T x 4  matrix of indices of tetrahedral corners
  %   type  string containing type of mass matrix to compute
  %     'barycentric': diagonal lumped mass matrix obtained by summing 1/3
  %       of volumes of surrounding tets
  %     Not yet supported:
  %       'full': full mass matrix for p.w. linear fem
  %       'voronoi'
  % Output:
  %   M  #V x #V matrix of cot weights 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: massmatrix
  %
  warning('Deprecated. Call massmatrix directly');
  M = massmatrix(V,T,type);
end

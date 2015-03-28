function [D] = div3(V,T)
  % DIV3 Compute the discrete divergence operator
  %
  % D = div(V,T)
  %
  % Inputs:
  %   V  #vertices by dim list of mesh vertex positions
  %   T  #elements by 3 list of mesh tet indices
  % Outputs:
  %   D  #vertices by #elements*3 divergence operator
  %

  warning('Deprecated. Call grad directly.');
  D = div(V,T);

  %G = grad3(V,T); 
  %vol = volume(V,T);
  %% Off by a factor 2?
  %D = -G'*repdiag(diag(sparse(vol)),3);
end


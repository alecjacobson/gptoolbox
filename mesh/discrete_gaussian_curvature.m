function k = discrete_gaussian_curvature(V,F)
  % DISCRETE_GAUSSIAN_CURVATURE Compute discrete gaussian curvature according
  % to (9) in "Discrete Differential-Geometry Operators for Triangulated
  % 2-Manifolds" [Meyer et al. 02] but without the inverse area term.
  % 
  % k = discrete_gaussian_curvature(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indies
  % Outputs:
  %   k  #V by 1 list of discrete gaussian curvature values
  %

  %K_G(x_i) = (2π - ∑θj)
  k = 2*pi - sparse(F,1,internalangles(V,F),size(V,1),1);



end

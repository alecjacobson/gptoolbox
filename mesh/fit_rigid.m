function [R,t] = fit_rigid(V,U)
  % FIT_RIGID fit a rigid transformation that best matches V to U in the least
  % squares sense:
  %   min ∑ || vi - R(ui+t) ||²
  %
  % Inputs:
  %   V  #V by dim list of points
  %   U  #V by dim list of points
  % Outputs:
  %   R  dim by dim rotation matrix transposed
  %   t  dim translation vector
  %
  % Example:
  %   [R,t] = fit_rigid(V,U);
  %   UU = bsxfun(@plus,U*R,t);

  % sum up outerproducts of vector from each vertex to center of mass in U and
  % V
  S = bsxfun(@minus,U,mean(U))'*bsxfun(@minus,V,mean(V));
  % Find closest rotation
  R = fit_rotation(S');

  % translation is just difference of center of mass (assuming uniform mass
  % distribution)
  t = mean(V)-mean(U)*R;

end

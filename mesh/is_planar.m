function [f,s] = is_planar(V,epsilon)
  % IS_PLANAR Return whether the 3D mesh is planar
  % 
  % [f,s] = is_planar(V)
  % [f,s] = is_planar(V,epsilon)
  %
  % Inputs:
  %   V  #V x 3 matrix of vertex coordinates
  %   epsilon  optional parameter used as threshold for smallest eigen value of
  %     pca, {1e-10}
  % Outputs:
  %   f  flag whether mesh is planar
  %   s  smallest principle component variance
  %
  %


  if ~exist('epsilon','var')
    epsilon = 1e-10;
  end


  if size(V,2) == 2 
    f = true;
    s = 0;
    return;
  end

  % This seems to have changed during matlab. For large V this is
  % infeasible.
  [c,l] = pcacov(V);

  s = min(l);


  f = s < epsilon;

end

function k = discrete_gaussian_curvature(V,F,varargin)
  % DISCRETE_GAUSSIAN_CURVATURE Compute discrete gaussian curvature according
  % to (9) in "Discrete Differential-Geometry Operators for Triangulated
  % 2-Manifolds" [Meyer et al. 02] but without the inverse area term.
  % 
  % k = discrete_gaussian_curvature(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indies
  % Optional:
  %     'ModifyBoundary' followed by true of {false} whether to subtract π
  %     at boundary vertices
  % Outputs:
  %   k  #V by 1 list of discrete gaussian curvature values
  %
  
  % Default value
  boundary = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'ModifyBoundary'}, {'boundary'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  %K_G(x_i) = (2π - ∑θj)
  k = 2*pi - sparse(F,1,internalangles(V,F),size(V,1),1);
  
  % Boundary vertices
  %K_G(x_i) = (π - ∑θj)
  if boundary
    b = outline(F);
    b = unique(b(:));
    k(b) = k(b) - pi;
  end

end

function [CV,CF,OV,OF] = closing(V,F,sigma,varargin)
  

  epsilon = [];
  st = 'pseudonormal';
  ct = 'marching_cubes';
  gs = 100;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Epsilon','SignedDistanceType','ContouringMethod','GridSize'}, ...
    {'epsilon','st','ct','gs'});
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
  assert(strcmp(ct,'marching_cubes'),'alt contouring not really supported yet');

  if isempty(epsilon)
    h = max(max(V)-min(V))/gs;
    epsilon = sign(sigma)*max(0.01*sigma,0.1*h);
  end

  [OV,OF] = signed_distance_isosurface( ...
    V,F, ...
    'Level',(sigma+epsilon), ...
    'SignedDistanceType',st,'ContouringMethod','marching_cubes','GridSize',gs);
  if isempty(OF)
    CV = [];
    CF = [];
  else
    [CV,CF] = signed_distance_isosurface( ...
      OV,OF, ...
      'Level',-sigma, ...
      'SignedDistanceType',st,'ContouringMethod','marching_cubes','GridSize',gs);
  end


end

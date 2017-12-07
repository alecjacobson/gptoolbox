function [V,F] = subdivided_sphere(iters,varargin)
  % SUBDIVIDED_SPHERE Generate a sphere by iteratively subdividing a
  % icosahedron in-plane and normalizing vertex locations to lie on the sphere.
  %
  % [V,F] = subdivided_sphere(iters,varargin)
  %
  % Inputs:
  %   iters  number of subdivision iterations (0 produces 12-vertex
  %     icosahedron)
  %   Optional:
  %     'Radius'  followed by scalar radius (multiplied against final V) {1}
  %     'SubdivisionMethod' followed by either:
  %       'loop'
  %       'sqrt3'
  %       {'upsample'}
  % Outputs:
  %   V  #V by 3 list of mesh vertices
  %   F  #F by 3 list of face indices into V
  % 
  % See also: upsample, loop
  % 

  % default values
  radius = 1;
  subdivision_method = 'upsample';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Radius','SubdivisionMethod'},{'radius','subdivision_method'});
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

  % Compute the 12 vertices
  phi = (1+sqrt(5))/2;  % Golden ratio
  V = [0   1  phi 
          0  -1  phi 
          0   1 -phi 
          0  -1 -phi 
          1  phi  0  
        -1  phi  0  
          1 -phi  0  
        -1 -phi  0  
          phi 0   1  
        -phi 0   1  
          phi 0  -1  
        -phi 0  -1];
  % Scale to required radius
  V = V/(sqrt(1+phi^2));
  % Define the adjacency matrix
  F = [1  2  9
        1  9  5
        1  5  6
        1  6  10
        1  10 2
        2  7  9
        9  7  11
        9  11 5
        5  11 3
        5  3  6
        6  3  12
        6  12 10
        10 12 8
        10 8  2
        2  8  7
        4  7  8
        4  8  12
        4  12 3
        4  3  11
        4  11 7];
  V = normalizerow(V);
  for iter = 1:iters
    switch subdivision_method
    case 'upsample'
      [V,F] = upsample(V,F);
    case 'loop'
      [V,F] = loop(V,F);
    case 'sqrt3'
      [V,F] = sqrt3(V,F);
    end
    V = normalizerow(V);

  end
  V = V*radius;

end


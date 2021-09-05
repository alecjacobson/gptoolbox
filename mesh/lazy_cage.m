function [dV,dF,b] = lazy_cage(V,F,m,varargin)
  % LAZY_CAGE  Build a cage that strictly encloses an input triangle soup.
  %
  % [dV,dF,b,IV,IF] = lazy_cage(V,F,m,varargin)
  % 
  % Inputs:
  %   V  #V by 3 list of input mesh vertex positions
  %   F  #F by 3 list of input mesh triangle indices into rows of V
  %   m  target number of faces of cage
  %   Optional:
  %     'DecimationMethod'  followed by 'Method' parameter to decimate_libigl
  %     'GridSize'  followed by size of grid to use for marching cubes step
  %       {ceil(sqrt(m)*3.17)}
  %     'MaxIter'  followed by maximum number of iterations {10}
  %     'MaxOffset'  followed by maximum offset value. Default value
  %       {norm(max(V)-min(V))/2} is _very_ conservative, but due to binary
  %       search usually just means 2-3 extra iterations.
  %     'Solid'  followed by whether the input should be treated as "solid",
  %       i.e. keep just the outer hull of each connected component or retain
  %       boundaries of internal cavities (perhaps usefull for enclosing
  %       shells). Solids with legitimate internal cavities will not currently
  %       be handled well. {true}
  % Outputs:
  %   dV  #dV by 3 list of output vertex positions
  %   dF  #dF by 3 list of output triangle indices into rows of dV
  %   b  smallest feasible offset found
  %

  % Larger → slower and more accurate
  gs = ceil(sqrt(m)*3.17);
  init_max_b = norm(max(V)-min(V))/2;
  decimation_method = 'qslim';
  max_iter = 10;
  solid = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'GridSize','DecimationMethod','MaxOffset','MaxIter','Solid'}, ...
    {'gs','decimation_method','init_max_b','max_iter','solid'});
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

  % binary search on offset parameter
  bounds = [0 init_max_b];
  for iter = 1:max_iter
    progressbar(iter,max_iter);
    b = mean(bounds);
    % Mesh offset
    [IV,IF] = signed_distance_isosurface(V,F,'SignedDistanceType','unsigned','Level',b,'GridSize',gs);
    % only keep outer shells
    if solid
      [~,C] = connected_components(IF);
      [~,vols] = arrayfun(@(c) centroid(IV,IF(C==c,:)),(1:max(C))','UniformOutput',0);
      [IV,~,~,IF] = remove_unreferenced(IV,IF(cell2mat(vols(C))>0,:));
    end
    if ~isempty(intersect_other(IV,IF,V,F,'FirstOnly',true))
      %fprintf('IV,IF intersects V,F\n');
      bounds(1) = b;
      continue;
    end
    [dV,dF,dJ] = decimate_libigl(IV,IF,m,'Method','naive');
    if ~isempty(intersect_other(dV,dF,V,F,'FirstOnly',true))
      %fprintf('dV,dF intersects V,F\n');
      bounds(1) = b;
      continue;
    end
    % self-union to handle any new self-intersections
    [dV,dF] = mesh_boolean(dV,dF,[],[],'union');
    % success
    bounds(2) = b;
  end
end

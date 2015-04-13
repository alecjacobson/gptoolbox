function N = per_vertex_normals(V,F,varargin)
  % PER_VERTEX_NORMALS  Compute per-vertex (area-weighted) normals over a mesh
  % (V,F)
  %
  % N = per_vertex_normals(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     'Weighting' followed by one of the following:
  %        {'area'} weighting
  %        'uniform' weighting
  %        'angle' weighting
  % Outputs:
  %   N  #V by 3 list of vertex normals
  %

  weighting = 'area';
  U = [];
  R = [];
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Weighting'},{'weighting'});
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

  FN = normalizerow(normals(V,F)+eps);
  switch weighting
  case 'uniform'
    W = ones(size(F,1),3);
  case 'angle'
    W = internalangles(V,F);
  case 'area'
    W = repmat(doublearea(V,F),1,3);
  end
  W = sparse(F(:),repmat(1:size(F,1),1,3),W,size(V,1),size(F,1));
  N = normalizerow(W*FN);

end

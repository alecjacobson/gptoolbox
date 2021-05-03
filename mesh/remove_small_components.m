function [U,G,I,J] = remove_small_components(V,F,varargin)
  % REMOVE_SMALL_COMPONENTS
  %
  % [U,G,I,J] = remove_small_components(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into rows of V
  %   Optional:
  %     'Tol' followed by minimum volume (or surface area) to use OR the word
  %       'max' to keep only biggest component {0.0001*total_vol}
  %     'Volumetric' followed by whether to use volume or surface area
  % Outputs:
  %   U  #U by 3 list of vertex positions
  %   G  #G by 3 list of face indices into rows of U
  %   I  #V by 1 list of indices such that: G = I(F)
  %   J  #G by 1 list of indices into F
  %

  tol = [];
  volumetric = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Tol','Volumetric'}, {'tol','volumetric'});
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


  if isempty(tol)
    if volumetric
      [~,total] = centroid(V,F);
    else
      total = sum(doublearea(V,F))*0.5;
    end
    tol = 0.0001*total;
  end


  [~,C] = connected_components(F);
  nc = max(C);
  val = zeros(nc,1);
  for i = 1:nc
    Fi = F(C==i,:);
    if volumetric
      [~,val(i)] = centroid(V,Fi);
    else
      val(i) = 0.5*sum(doublearea(V,Fi));
    end
  end

  if isnumeric(tol)
    J = find(ismember(C,find(val>tol)));
  else 
    assert(strcmp(tol,'max'));
    [~,max_c] = max(val);
    J = find(C==max_c);
  end
  F = F(J,:);
  [U,I] = remove_unreferenced(V,F);
  G = I(F);
end

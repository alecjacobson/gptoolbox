function [U,G] = shell(V,F,th,varargin)
  % SHELL Compute a thin shell around a triangle mesh (V,F) with thickness th.
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   th  desired thickness of shell
  %   Optional:
  %     'Normals' followed by #V by 3 list of vertex normals to use to
  %     construct shell (will be multiplied by th): {area-weighted}
  % Outputs:
  %   U  2*#V by 3 list of output mesh vertex positions: V always comes first
  %   G  #G by 3 list of triangle indices into U
  %

  N = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Normals'}, ...
    {'N'});
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

  if isempty(N)
    N = normalizerow(per_vertex_normals(V,F));
  end

  % scale by thickness
  N = N*th;

  O = outline(F);
  n = size(V,1);
  U = [V;V+N];
  G = [ ...
    fliplr(F); ...
    n+F; ...
    0+O(:,1) n+O(:,[2 1]); ...
    n+O(:,2) 0+O(:,1:2)];

end

function [U,data] = fast_mass_springs(V,E,b,bc,varargin)
  % FAST_MASS_SPRINGS Solve for the mass springs deformation subject to hard
  % constraints, following "Fast Simulation of Mass-Spring Systems" [Liu et al.
  % 2013].
  %
  % U = fast_mass_springs(V,E,b,bc)
  %
  % Inputs:
  %   V  #V by dim list curve vertex positions
  %   E  #E by 2 list of *edge* indices or other simplices and edges will be
  %     derived.
  %   b  #b list of indices into V of fixed vertices
  %   bc  #b by dim list of positions for fixed vertices
  %    Optional:
  %       'EdgeLengths' followed by #E list of rest edge lengths
  %       'Data' computed min_quad_with_fixed data
  %       'U0'  followed by #V by dim initial guess
  % Outputs:
  %   U  #V by dim list of deformed curve vertex positions
  %   data  struct output of min_quad_with_fixed
  %

  if size(E,2)>2
    E = edges(E);
  end

  % default values
  R = [];
  data = [];
  % Initial guess 
  U = V;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'EdgeLengths','Data','U0'}, ...
    {'R','data','U'});
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

  % number of vertices
  m = size(V,1);
  % number of springs
  s = size(E,1);

  if isempty(data)
    % compute rest edge norms
    if isempty(R)
      data.R = sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
    end
    % Inversely proportional to length (then squared)
    K = data.R.^-0.5;
    %% Uniform for now
    %K = ones(s,1);
    % Spring signed incidence matrix
    A = sparse(E,[1:s;1:s]',[K,-K],m,s);
    % "Laplacian"
    data.L = A*A';
    S = diag(sparse(K));
    data.J = A*S;
    data.mqwf = [];
  end

  iter = 0;
  max_iter = 100;
  while true
    % Fix U and find D ("local" step)
    % recompute edge vectors
    D = U(E(:,1),:)-U(E(:,2),:);
    % normalize
    D = bsxfun(@times,D,data.R./sqrt(sum(D.^2,2)));
    % Fix D and find U ("global" step)
    [U,data.mqwf] = min_quad_with_fixed(0.5*data.L,-data.J*D,b,bc,[],[],data.mqwf);
    iter = iter + 1;
    if iter == max_iter
      break;
    end
  end
end

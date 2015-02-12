function U = fast_mass_springs(V,E,b,bc,varargin)
  % FAST_MASS_SPRINGS Solve for the mass springs deformation subject to hard
  % constraints, following "Fast Simulation of Mass-Spring Systems" [Liu et al.
  % 2013].
  %
  % U = fast_mass_springs(V,E,b,bc)
  %
  % Inputs:
  %   V  #V by 2 list curve vertex positions
  %   E  #E by 2 list of *edge* indices or other simplices and edges will be
  %     derived.
  %   b  #b list of indices into V of fixed vertices
  %   bc  #b by 2 list of positions for fixed vertices
  %    Optional:
  %       'EdgeLengths' followed by #E list of rest edge lengths
  % Outputs:
  %   U  #V by 2 list of deformed curve vertex positions
  %

  if size(E,2)>2
    E = edges(E);
  end

  % default values
  R = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'EdgeLengths'}, ...
    {'R'});
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

  % number of dimensions
  dim = size(V,2);
  % number of vertices
  m = size(V,1);
  % number of springs
  s = size(E,1);
  % Initial guess 
  U = V;
  % compute rest edge norms
  if isempty(R)
    R = sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
  end
  % Inversely proportional to length (then squared)
  K = R.^-0.5;
  %% Uniform for now
  %K = ones(s,1);
  % Spring incidence matrix
  A = sparse(E,[1:s;1:s]',[K,-K],m,s);
  % dim by dim identity
  I = speye(dim,dim);
  % "Laplacian"
  L = repdiag(A*A',dim);
  S = diag(K);
  J = repdiag(A*S,dim);
  % repeat boundary for dimensions
  b_dim = [b(:);m+b(:)];

  iter = 0;
  max_iter = 10;
  mqwf = [];
  while true
    % Fix U and find D ("local" step)
    % recompute edge vectors
    D = U(E(:,1),:)-U(E(:,2),:);
    % normalize
    D = bsxfun(@times,D,R./sqrt(sum(D.^2,2)));
    % Fix D and find U ("global" step)
    [U,mqwf] = min_quad_with_fixed(0.5*L,-J*D(:),b_dim,bc(:),[],[],mqwf);
    U = reshape(U,size(V));
    
    iter = iter + 1;
    if iter == max_iter
      break;
    end
  end
end

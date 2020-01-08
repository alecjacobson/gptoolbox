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
  % quadratic term
  Q = sparse(0);
  % linear term
  l = sparse(0);
  % stiffness
  k = 1;
  % Initial guess 
  U = V;
  max_iter = 100;
  Aeq = [];
  Beq = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Aeq','Beq','EdgeLengths','Stiffness','Data','Q','l','MaxIter','U0'}, ...
    {'Aeq','Beq','R','k','data','Q','l','max_iter','U'});
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
  assert((numel(k) == 1)||(numel(k) == s));

  if isempty(data)
    % compute rest edge norms
    if isempty(R)
      data.R = sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
    else
      data.R = R;
    end
    % Inversely proportional to length (then squared)
    %
    % Huh? In the paper, A is really just the incidence matrix with ±1
    K = sqrt(k).*data.R.^-0.5;
    %% Uniform for now
    %K = ones(s,1);
    % Spring signed incidence matrix
    A = sparse(E,[1:s;1:s]',[K,-K],m,s);
    % "Laplacian"
    data.L = A*A';
    S = diag(sparse(K));
    % Suspicious that sqrt(k) shows up here but k shows up in L
    data.J = A*S;
    if iscell(Aeq)
      data.mqwf = cell(size(Aeq));
    else
      data.mqwf = [];
    end
  end

  iter = 0;
  while true
    % Fix U and find D ("local" step)
    % recompute edge vectors
    D = U(E(:,1),:)-U(E(:,2),:);
    % normalize
    D = bsxfun(@times,D,data.R./sqrt(sum(D.^2,2)));
    % Fix D and find U ("global" step)
    if iscell(Aeq)
      U = zeros(size(V));
      for c = 1:size(V,2)
        [U(:,c),data.mqwf{c}] =  min_quad_with_fixed( ...
          0.5*data.L+Q,-data.J*D(:,c)+l(:,c),b,bc(:,c),Aeq{c},Beq{c},data.mqwf{c});
      end
    else
      [U,data.mqwf] = min_quad_with_fixed(0.5*data.L+Q,-data.J*D+l,b,bc,Aeq,Beq,data.mqwf);
    end

    iter = iter + 1;
    if iter == max_iter
      break;
    end
  end
end

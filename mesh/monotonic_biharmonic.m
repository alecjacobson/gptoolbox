function [W,f] = monotonic_biharmonic(varargin)
  % MONOTONIC_BIHARMONIC  Compute "monotonic biharmonic" coordinates, using
  % quadratic programming optimization. Each weight funciton minimizes
  % laplacian energy subject to monotonicity constraints taken from harmonic
  % coordinates
  %
  % W = monotonic_biharmonic(V,F,b,bc,'ParameterName',ParameterValue)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %    Optional:
  %      'k'  followed by integer, solve k-harmonic problem {2}
  %      'QuadProgParam'  followed by mosek param struct, see optimset. If not
  %        given then *some* default options are used
  %      'MonotonicityGraph' followed A a #handles cell of #n by #n 
  %        monotonicity graph (sparse adjacency matrix) with A{k}(i,j) !=
  %        0 if we want that W(i,k) <= W(j,k), default is to use graph
  %        from harmonic coordinates
  %      'MonotonicityTGraph' followed AT a #handles cell of #F by #F set of
  %        constraints on triangle (averages)
  %      'Data' followed by #W by #handles data term to be
  %        'as-close-as-possible', {[]}
  %      'DataWeight' follow by weight of Data energy term, {1}
  %      'Algorithm' followed by solver name {'interior-point'}:
  %         'interior-point' 
  %         'active-set' Alec's active set solver
  %           (min_quad_with_fixed_active_set)
  %      'OptType'
  %         optimize using {'conic'} or 'quad' or 'legacy-quad'
  %      'SubGraph' followed by true or false, whether to iteratively resolve
  %        using exponentially larger subgraphs of the monotonicity graph,
  %        {false}
  %      'IterateConstraints' Iterate on monotonicity constraints
  %      'ExtremaBC'  Followed by bc to be used for computing monotonicitygraph
  %      0 means min, 1 means max
  %      'MaxIter' maximum number of iterations {inf}
  %      'Quiet'  followed by 'echo(0)' or '' to be quite or loud
  %      'MonotonicityEpsilon' followed by epsilon used for enforcing
  %        monotonicity. {0} Increase this value if constraints are being
  %        violated.
  %      'LEQScale' followed by scale applied to linear inequality constraints
  %        to help mosek know that our constraints are not to be violated {1000}.
  %        This is probably tied to the implicit scaling of the data and 
  %        smoothness terms. Increase this value if constraints are being
  %        violated
  % Output:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions, biharmonic_bounded
  %

  V = varargin{1};
  F = varargin{2};
  b = varargin{3};
  bc = varargin{4};

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);

  % default options
  k = 2;
  algorithm = 'interior-point';
  opt_type = 'conic';
  subgraph = false;
  iterate_constraints = false;
  extrema_bc = bc;
  max_iter = inf;
  % check for mosek and set its parameters
  [param,mosek_exists] = default_mosek_param();
  if mosek_exists
    opt_type = 'conic';
  else
    opt_type = 'legacy-quad';
  end
  %param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 0;
  %param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 0;

  S = [];
  data_weight = 1;
  A = {};
  AT = cell(m,1);
  quiet = 'echo(0)';
  monotonicity_epsilon = 0;
  leq_scale = 1;

  % parse optional inputs
  ii = 5;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Quiet'
      ii = ii + 1;
      assert(ii<=nargin);
      quiet = varargin{ii};
    case 'MonotonicityGraph'
      ii = ii + 1;
      assert(ii<=nargin);
      A = varargin{ii};
    case 'MonotonicityTGraph'
      ii = ii + 1;
      assert(ii<=nargin);
      AT = varargin{ii};
    case 'Data'
      ii = ii + 1;
      assert(ii<=nargin);
      S = varargin{ii};
    case 'DataWeight'
      ii = ii + 1;
      assert(ii<=nargin);
      data_weight = varargin{ii};
    case 'k'
      ii = ii + 1;
      assert(ii<=nargin);
      k = varargin{ii};
    case 'QuadProgParam'
      ii = ii + 1;
      assert(ii<=nargin);
      param = varargin{ii};
    case 'Algorithm'
      ii = ii + 1;
      assert(ii<=nargin);
      algorithm = varargin{ii};
    case 'OptType'
      ii = ii + 1;
      assert(ii<=nargin);
      opt_type = varargin{ii};
    case 'SubGraph'
      ii = ii + 1;
      assert(ii<=nargin);
      subgraph = varargin{ii};
    case 'IterateConstraints'
      ii = ii + 1;
      assert(ii<=nargin);
      iterate_constraints = varargin{ii};
    case 'ExtremaBC'
      ii = ii + 1;
      assert(ii<=nargin);
      extrema_bc = varargin{ii};
    case 'MaxIter'
      ii = ii + 1;
      assert(ii<=nargin);
      max_iter = varargin{ii};
    case 'MonotonicityEpsilon'
      ii = ii + 1;
      assert(ii<=nargin);
      monotonicity_epsilon = varargin{ii};
    case 'LEQScale'
      ii = ii + 1;
      assert(ii<=nargin);
      leq_scale = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii+1;
  end

  % Build discrete laplacian and mass matrices used by all handles' solves
  if(size(F,2)==4)
    fprintf('Solving over volume...\n');
    L = cotmatrix3(V,F);
    M = massmatrix3(V,F,'barycentric');
    % Build pieces of sqrt(L) for odd powers
    if mod(k,2) ~= 0
      % area is really volume
      dblA = abs(volume(V,F));
      TA = repdiag(diag(sparse(dblA)),size(V,2));
      G = grad3(V,F);
    end
  else
    L = cotmatrix_embedded(V,F);
    M = massmatrix_embedded(V,F,'voronoi');
    % Build pieces of sqrt(L) for odd powers
    if mod(k,2) ~= 0
      dblA = doublearea(V,F);
      TA = repdiag(diag(sparse(dblA)/2),size(V,2));
      G = grad(V,F);
    end
  end
  % Normalize mass matrix (important!)
  M = M./max(abs(diag(M)));

  % build quadratic coefficient matrix (bilaplacian operator) 
  switch opt_type
  case 'legacy-quad'
    Q = -L;
    % Construct quadratic coefficients matrix
    for ii = 2:k
      Q = -L*(M\Q);
    end
    % Symmetrize for high-order
    if k >= 3
      Q = (Q+Q')*0.5;
    end
    save('mbQ.mat','Q');
  case {'quad','conic'}
    % build possibly (odd k) rectangular sqrt of quadratic coefficient matrix
    if mod(k,2)==0
      FF = -L;
      FFstr = '-L';
      % even k
      for ii = 2:(k/2)
        FF = -L*(M\FF);
        FFstr = ['-L(M\' FFstr ')'];
      end
      FF = sqrt(M)\FF;
      FFstr = ['sqrt(M)\(' FFstr ')'];
    else
      % odd k
      FF = speye(n,n);
      FFstr = '';
      for ii = 1:((k-1)/2)
        FF = M\(-L*FF);
        FFstr = ['M\(-L ' FFstr ')'];
      end
      FF = sqrt(TA)*G*FF;
      FFstr = ['sqrt(TA)G(' FFstr ')'];
    end
  otherwise
    error(['Unsupported OptType: ' opt_type]);
  end

  % No linear terms
  l = zeros(n,m);

  if ~isempty(S)
    % build quadratic and linear data term
    if numel(data_weight) == 1
      data_weight = repmat(data_weight,n,1);
    end
    % L2 norm scaling
    data_scale = 1/sqrt(sum(S(:).^2));
    switch opt_type
    case 'legacy-quad'
      % Mosek assumes quadratic terms will be multiplied by 1/2
      Q = Q+2*data_scale*diag(sparse(data_weight));
    case {'quad','conic'}
      FF = [FF;diag(sqrt(2)*sparse(sqrt(data_scale*data_weight)))];
    end
    l = bsxfun(@plus,l,-2*bsxfun(@times,data_scale*data_weight,S));
  end

  if isempty(A)
    % initial monotonic guess 
    assert(all(extrema_bc(:)==0 | extrema_bc(:)==1));
    WH = kharmonic(V,F,b,extrema_bc);
    A = cell(m,1);
    % loop over handles
    for i = 1:m
      A{i} = monotonicity_matrix(WH(:,i),F,'OneInOneOut',true,'EnforceExtrema',true);
    end
  elseif isnumeric(A)
    assert(~iterate_constraints);
    A = {A};
  end
  if isnumeric(AT)
    AT = {AT};
  end

  tic;
  % loop over handles
  for i = 1:m
    iter = 1;
    fullgraph = false;
    while(iter <= max_iter)
      Ai = A{i};
      ATi = AT{i};
      % Unpublished idea to only enforce some subset of the constraints, you
      % may safely ignore this section
      if subgraph
        assert(isempty(ATi));
        [AiI,AiJ,AiV] = find(Ai);
        step_base = 2;
        num_edges = round(step_base^(1.1*(iter+1))-1);
        if num_edges < numel(AiI)
          fprintf('Num edges in subgraph: %d\n',num_edges);
          randorder = randperm(numel(AiI));
          sel = randorder(1:num_edges);
          Ai = sparse(AiI(sel),AiJ(sel),AiV(sel),size(Ai,1),size(Ai,2));
        else
          num_edges  = numel(AiI);
          fprintf('Num edges in fullgraph: %d\n',numel(AiI));
          fullgraph = true;
        end
      end
      % monotonicity constraints, as less than or equals constraints
      [Aleq,ATleq] = inequality_constraints_from_graphs(F,Ai,ATi);
      Aleq = [Aleq;ATleq];
      nleq = size(Aleq,1);
      bleq = zeros(nleq,1);

      % boundary conditions
      Aeq = speye(n,n);
      Aeq = Aeq(b,:);

      % Optimize according to algorithm and opt_type choice
      switch algorithm
      case 'interior-point'
        % Constant bounds seem to help performance
        lx = min(bc(:,i))*ones(n,1);
        ux = max(bc(:,i))*ones(n,1);
        switch opt_type
        case 'legacy-quad'
          if(mosek_exists)
            fprintf('Quadratic optimization using mosek...\n');
          else
            fprintf('Quadratic optimization using matlab...\n');
          end
          fprintf( [ ...
            '  minimize:     x''LM\\Lx\n' ...
            'subject to: xi <= xj, ???j, j parent to i\n' ]);
          [x,fval,err] = quadprog(Q,l(:,i),leq_scale*Aleq,-monotonicity_epsilon+leq_scale*bleq, ...
            Aeq,bc(:,i),lx,ux,[],param);
          if(err ~= 1)
            fprintf([...
              '----------------------------------------------------------\n' ...
              'ERROR ('  num2str(err) ',' num2str(fval) '):' ...
              ' solution may be inaccurate...\n' ...
              '----------------------------------------------------------\n' ...
              ]);
          end
        case {'quad','conic'}
          fprintf( [ ...
            '  minimize:     ||Fx||^2\n' ...
            'subject to: xi <= xj, ???j, j parent to i\n' ]);
          [x,f] = conic( ...
            FF,l(:,i), ...
            [],[], ...
            leq_scale*Aleq,-monotonicity_epsilon+leq_scale*bleq, ...
            lx,ux, ...
            b,bc(:,i), ...
            'Param',param, ...
            'Quiet',quiet, ...
            'OptType',opt_type);
        end
        % number of violated constraints
        diffleq = Aleq*x - bleq;
        vc = sum(diffleq>=0);
        if vc>0
          warning(['%d violated constraints, ' ...
            'try increasing MonotonicityEpsilon by max violation: +%g ' ...
            'or increasing LEQScale '], ...
            vc,max(diffleq(diffleq>0)));
        end
      % Unpublished home-brewed active set solver, likely this will not work.
      % Safe to ignore
      case 'active-set'
        fprintf('Quadratic optimization using active set...\n');
        [x,preF] = min_quad_with_fixed_active_set(Q,zeros(n,1),b,bc(:,i),[],[],[],[]);
      otherwise
        error(['Unsupported algorithm type: ''' algorithm '''']);
      end
      % set weights to solution in weight matrix
      W(:,i) = x(1:n);
      fprintf('Lap time: %gs\n',toc);

      if ~subgraph && ~iterate_constraints
        break;
      elseif subgraph
        % Again, unpublished and safe to ignore this section
        if size(V,2) == 2
          Z = W(:,i);
        else
          Z = V(:,3);
        end
        t = trisurf(F,V(:,1),V(:,2),Z,'FaceLighting','phong','FaceColor','interp','CData',W(:,i),'EdgeColor','none');
        colormap(jet(20));
        mes = mark_extrema(t,W(:,i));
        drawnow;
        if max_iter == inf
          %input(['iter: ' num2str(iter)]);
        end
      elseif iterate_constraints
        % Rebuild constraints based on previous solution
        A{i} = monotonicity_matrix(W(:,i),F,'OneInOneOut',true,'EnforceExtrema',true);
        if max_iter == inf
          %input(['iter: ' num2str(iter)]);
        end
      end

      % Unpublished, safe to ignore
      if fullgraph
        break;
      end 

      % Unpublished, safe to ignore
      if iter == max_iter
        fprintf('Reached maximum iterations (%d) with convergence\n',max_iter);
        break;
      end
      iter = iter + 1;
    end
  end
  t = toc;
  % Print some time information
  fprintf('Total elapsed time: %gs\n',t);
  fprintf('Average time per handle: %gs\n',t/m);

end

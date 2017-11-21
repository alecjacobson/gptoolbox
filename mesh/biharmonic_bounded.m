function W = biharmonic_bounded(V,F,b,bc,varargin)
  % BIHARMONIC_BOUNDED Compute biharmonic bounded coordinates, using quadratic
  % optimizer
  %
  % W = biharmonic_bounded(V,F,b,bc)
  % W = biharmonic_bounded(V,F,b,bc,'ParameterName',ParameterValue)
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
  %      'OptType'  of optimizer to use {best available}:
  %        'quad'
  %        'least-squares'
  %        'conic'
  %        'active-set'
  %      'POU'  true or false, enforce partition of unity explicitly {false}
  %      'Low'  lower bound {0}
  %      'Up'  upper bound {1}
  %      'ShapePreserving'  #V rigidity mask, where 0 means not enforced, 1 means strongly
  %        enforced
  %  
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions
  %

  function COUNT = verbose_fprintf(varargin)
    if verbose
      COUNT = fprintf(varargin{:});
    end
  end

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);


  % default options
  pou = false;
  k = 2;
  % check for mosek and set its parameters
  [param,mosek_exists] = default_quadprog_param();
  param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-12;
  if mosek_exists
    opt_type = 'conic';
  else
    opt_type = 'quad';
  end
  % default bounds
  low = 0;
  up = 1;
  R = sparse(n,1);



  % default values
  verbose = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Verbose','W0','k','QuadProgParam','Low','Up','POU','OptType','ShapePreserving'}, ...
    {'verbose','W0','k','param','low','up','pou','opt_type','R'});
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

  assert(k <= 2 || strcmp(opt_type,'quad'), ...
    'For large k, OptType should be quad.');

  % Build discrete laplacian and mass matrices used by all handles' solves
  L = cotmatrix(V,F);
  M = massmatrix(V,F);
  % NORMALIZE MASSMATRIX (THIS IS IMPORTANT!!)
  M = M./max(abs(diag(M)));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP SOLVER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if any(R(:)) && ~strcmp(opt_type,'quad')
    error( [ ...
      'Enforcing shape preserving regions only supported in' ...
      ' conjunction with opt_type=''quad''']);
  end

  %% Shape preserving laplacian
  %Lsp = bsxfun(@times,L,sparse(R)*shape_preserving_weight);
  %spy(Lsp)
  shape_preserving_weight = 100;
  E = edges(F);
  ne = size(E,1);
  % Build gradient
  if any(R(:))
    G = sparse( ...
      [1:ne 1:ne]', ...
      [E(:,1);E(:,2)], ...
      shape_preserving_weight*[mean(R(E),2);-mean(R(E),2)], ...
      ne, ...
      n);
    Lsp = G' * G;
  else
    G = sparse(0,n);
    Lsp = sparse(n,n);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP PROBLEM AND SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(pou)
    % Enforce partition of unity as explicity constraints: solve for weights
    % of all handles simultaneously
    if(strcmp(opt_type,'quad'))
      % biharmonic system matrix
      Qi = -L;
      for ii = 2:k
        Qi = -L*(M\Qi);
      end
      if k>3
        Qi = (Qi + Qi')*0.5;
      end
      Qi = Qi + Lsp;
      Q = sparse(m*n,m*n);
      % Q is sparse matrix with Qi along diagonal
      for ii = 1:m
        d = (ii - 1)*n + 1;
        Q(d:(d + n-1), d:(d + n-1)) = Qi;
      end
      % linear constraints: partition of unity constraints and boundary
      % conditions
      PA = repmat(speye(n,n),1,m);
      Pb = ones(n,1);
      % boundary conditions
      BCAi = speye(n,n);
      BCAi = BCAi(b,:);
      BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
      % BCA is sparse matrix with BCAi along diagonal
      for ii = 1:m
        di = (ii - 1)*size(BCAi,1) + 1;
        dj = (ii - 1)*size(BCAi,2) + 1;
        BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
      end
      BCb = bc(:);
      % set bounds
      ux = up.*ones(m*n,1);
      lx = low.*ones(m*n,1);
      if(mosek_exists)
        verbose_fprintf('Quadratic optimization using mosek...\n');
      else
        verbose_fprintf('Quadratic optimization using matlab...\n');
      end
      verbose_fprintf( [ ...
        '  minimize:     x''LM\\Lx\n' ...
        'subject to: %g <= x <= %g, ???_i xi = 1\n'], ...
        low,up);
      tic;
      W = quadprog(Q,zeros(n*m,1),[],[],[PA;BCA],[Pb;BCb],lx,ux,[],param);
      t = toc;
      verbose_fprintf('Total elapsed time: %gs\n',t);
      W = reshape(W,n,m);
    else
      error( [ ...
        'Enforcing partition of unity only support in conjunction with ' ...
        'opt_type=''quad''']);
    end
  else
    % Drop partition of unity constraints, solve for weights of each handle
    % independently then normalize to enforce partition of unity
    switch opt_type
    case {'quad','active-set'}
      % build quadratic coefficient matrix (bilaplacian operator)
      Q = -L;
      for ii = 2:k
        Q = -L*(M\Q);
      end
      if k>=3
        Q = (Q + Q')*0.5;
      end
      Q = Q + Lsp;
      % set bounds
      ux = up.*ones(n,1);
      lx = low.*ones(n,1);
    case 'least-squares'
      % solve same problem but as least-squares problem see mosek documention
      % for details
      I = speye(n);
      Z = sparse(n,n);
      Q = [Z,Z;Z,I];
      assert(k == 2);
      F = [sqrt(M)\L;G];
      c = zeros(n,1);
      B = [F,-I];
      ux = [up.*ones(n,1) ;  Inf*ones(n,1)];
      lx = [low.*ones(n,1); -Inf*ones(n,1)];
    case 'conic'
      % solve same problem but as conic problem see mosek documention for
      % details
      %
      % Variable names are consistent with mosek doc 7.9.1:
      % [x z t]
      assert(k==2);
      F = sqrt(M)\L;
      prob.c = [zeros(2*n,1); 1];
      I = speye(n);
      prob.a = [F,-I,zeros(n,1)];
      prob.blc = zeros(n,1);
      prob.buc = zeros(n,1);
      prob.bux = [ up.*ones(n,1);  Inf*ones(n,1);  Inf];
      prob.blx = [ low.*ones(n,1); -Inf*ones(n,1); 0];
      prob.cones = cell(1,1);
      prob.cones{1}.type = 'MSK_CT_QUAD';
      t_index = 2*n +1;
      z_indices = (n+1):(2*n);
      prob.cones{1}.sub = [t_index z_indices];
    otherwise
      error('Bad opt_type');
    end

    % number of handles
    m = size(bc,2);
    % allocate space for weights
    W = zeros(n,m);
    tic;
    if strcmp(opt_type,'active-set')
      if ~exist('W0') || isempty(W0)
        verbose_fprintf('Initial guess for active set...\n');
        W = min_quad_with_fixed(Q,[],b,bc);
        verbose_fprintf('Lap time: %gs\n',toc);
      else
        W = W0;
      end
    end
    % loop over handles
    for i = 1:m
      switch opt_type
      case 'active-set'
        verbose_fprintf('Quadratic optimization using active set...\n');
        verbose_fprintf( [ ...
          '  minimize:     x''LM\\Lx\n' ...
          'subject to: %g <= x <= %g\n' ], ...
          low,up);
        AS = [];
        AS.Z0 = W(:,i);
        x = min_quad_with_fixed_active_set(Q,[],b,bc(:,i),[],[],[],[],lx,ux,AS);
      case 'quad'
        % enforce boundary conditions via lower and upper bounds
        %lx(b) = bc(:,i);
        %ux(b) = bc(:,i);
        Aeq = speye(n,n);
        Aeq = Aeq(b,:);
        if(mosek_exists)
          verbose_fprintf('Quadratic optimization using mosek...\n');
        else
          verbose_fprintf('Quadratic optimization using matlab...\n');
        end
        verbose_fprintf( [ ...
          '  minimize:     x''LM\\Lx\n' ...
          'subject to: %g <= x <= %g\n' ], ...
          low,up);
        % if mosek is not available, then matlab will complain that sparse
        % matrices are not yet supported...
        [x,fval,err] = quadprog(Q,zeros(n,1),[],[],Aeq,bc(:,i),lx,ux,[],param);
        if(err ~= 1)
          verbose_fprintf([...
            '----------------------------------------------------------\n' ...
            'ERROR ('  num2str(err) ',' num2str(fval) '):' ...
            ' solution may be inaccurate...\n' ...
            '----------------------------------------------------------\n' ...
            ]);
        end

      case 'least-squares'
        % enforce boundary conditions via lower and upper bounds
        lx(b) = bc(:,i);
        ux(b) = bc(:,i);
        verbose_fprintf('Quadratic optimization using mosek...\n');
        verbose_fprintf([ ...
          '  minimize:       z''z\n' ...
          '  subject to: M\\Lx - z = 0\n' ...
          '  and          %g <= x <= %g\n'], ...
          low,up);
        x = quadprog(Q,zeros(2*n,1),[],[],B,c,lx,ux,[],param);
      case 'conic'
        prob.bux(b) = bc(:,i);
        prob.blx(b) = bc(:,i);
        verbose_fprintf('Conic optimization using mosek...\n');
        verbose_fprintf([ ...
          '  minimize:         t\n' ...
          '  subject to: M\\Lx - z = 0,\n' ...
          '             t >= sqrt(z''z),\n' ...
          '               %f <= x <= %f\n'], ...
          low,up);
        [r,res]=mosekopt('minimize echo(0)',prob,param);
        % check for mosek error
        switch r
        case 4006
          warning(['MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg ...
            'The solution is probably OK, but ' ...
            'to make this error go away, increase: ' ...
            'MSK_DPAR_INTPNT_CO_TOL_REL_GAP' ...
            n]);
        case 10006
          warning(['MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg ...
            'The solution is probably OK, but ' ...
            'to make this error go away, decrease: ' ...
            'MSK_DPAR_INTPNT_CO_TOL_REL_GAP' ...
            n]);
        otherwise
          if r ~= 0
            error(['FATAL MOSEK ERROR. rcode: ' ...
              num2str(res.rcode) ' ' ...
              res.rcodestr ' ' ...
              res.rmsg]);
          end
        end
        % extract solution from result
        x = res.sol.itr.xx;
      end
      % set weights to solution in weight matrix
      W(:,i) = x(1:n);
      verbose_fprintf('Lap time: %gs\n',toc);
    end
    t = toc;
    verbose_fprintf('Total elapsed time: %gs\n',t);
    verbose_fprintf('Average time per handle: %gs\n',t/m);
  end

end

function [Z,f] = conic(varargin)
  % CONIC optimization wrapper for mosek, uses same prototype for conic and
  % quadprog 
  % minimize     ||Fx||² + c'x
  % subject to:  Aeq x = beq
  %              Aleq x ≤ bleq
  %              x ≤ lx
  %              x ≥ ux
  %              x(b) = bc
  % 
  % Z = conic(F,c,Aeq,beq,Aleq,bleq,lx,ux,b,bc,'ParameterName',ParameterValue)
  %
  % Inputs:
  %    F  #F by #x quadratic, least-squares matrix
  %    c  #x by 1 linear term
  %    Aeq  neq by #x  linear equality matrix
  %    beq  neq by 1  linear equality rhs
  %    Aleq  nleq by #x  linear inequality matrix
  %    bleq  neq by 1  linear inequality rhs
  %    lx  #x by 1 lower bounds
  %    ux  #x by 1 upper bounds
  %    b  #b list of boundary (fixed) indices in x
  %    bc  #b list of boundary (fixed) values for b
  %  Optional:
  %    'Quiet'  followed by 'echo(0)' or '' to be quite or loud
  %    'OptType' followed by 'conic' or 'quad'
  %    'Param' followed by mosek param struct
  % 

  FF = varargin{1};
  c = varargin{2};
  Aeq = varargin{3};
  beq = varargin{4};
  Aleq = varargin{5};
  bleq = varargin{6};
  lx = varargin{7};
  ux = varargin{8};
  b = varargin{9};
  bc = varargin{10};

  param = default_mosek_param();
  quiet = 'echo(0)';
  opt_type = 'conic';

  ii = 11;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Quiet'
      ii = ii + 1;
      assert(ii<=nargin);
      quiet = varargin{ii};
    case 'Param'
      ii = ii + 1;
      assert(ii<=nargin);
      param = varargin{ii};
    case 'OptType'
      ii = ii + 1;
      assert(ii<=nargin);
      opt_type = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end
  if isfield(param,'Diagnostics');
    param = rmfield(param,'Diagnostics');
  end

  % number of primary variables
  n = size(FF,2);
  nt = size(FF,1);
  % identity the size of rows in FF
  I = speye(nt,nt);
  % number of equality constraints
  neq = size(Aeq,1);
  % number of less than or equals constraints
  nleq = size(Aleq,1);

  % avoid warnings when Aeq and Aleq are empty
  if isempty(Aeq)
    Aeq = sparse(0,n);
  end
  if isempty(Aleq)
    Aleq = sparse(0,n);
  end

  switch opt_type
  case 'conic'
    % Follows: mosek 6.0's documentation "10.3.3.1. Quadratic objective and
    % constraints" but ignores quadratic constraints because we don't have any.

    % set up mosek problem
    prob = [];
    prob.c = [c;zeros(nt,1);1;0];
    prob.a = [ ...
      FF    -I             zeros(nt,2); ...
      Aeq  sparse(neq,nt)  zeros(neq,2)      ; ...
      Aleq sparse(nleq,nt) zeros(nleq,2)     ; ...
      ];
    prob.blc = [zeros(nt,1);beq;-Inf(nleq,1)];
    prob.buc = [zeros(nt,1);beq;bleq];
    prob.blx = [lx;-Inf(nt,1);0;1];
    prob.bux = [ux;Inf(nt,1);Inf;1];
    prob.cones = cell(1,1);
    prob.cones{1}.type = 'MSK_CT_RQUAD';
    v_index = n+nt+1;
    w_index = n+nt+2;
    t_indices = n+(1:nt);
    prob.cones{1}.sub = [v_index w_index t_indices];

    % Enforce fixed values
    prob.bux(b) = bc;
    prob.blx(b) = bc;

    fprintf('Conic optimization using mosek...\n');
  case 'conic-debug'
    prob = [];
    t_indices = n+(1:nt);
    prob.qosubi = t_indices;
    prob.qosubj = t_indices;
    prob.qoval = repmat(1,nt,1);
    prob.c = [c;zeros(nt,1)];
    prob.a = [ ...
      FF -I; ...
      Aeq  sparse(neq,nt)  ; ...
      Aleq sparse(nleq,nt) ; ...
      ];
    prob.blc = [zeros(nt,1);beq;-Inf(nleq,1)];
    prob.buc = [zeros(nt,1);beq;bleq];
    prob.blx = [lx;-Inf(nt,1)];
    prob.bux = [ux;Inf(nt,1)];

    % Enforce fixed values
    prob.bux(b) = bc;
    prob.blx(b) = bc;

    fprintf('Conic-debug optimization using mosek...\n');
  case 'quad'

    prob = [];
    [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(FF'*FF));
    prob.c = c;
    prob.a = [Aeq; Aleq];
    prob.blc = [beq;-Inf(nleq,1)];
    prob.buc = [beq;bleq];
    prob.blx = lx;
    prob.bux = ux;

    % Enforce fixed values
    prob.bux(b) = bc;
    prob.blx(b) = bc;
    fprintf('Quadratic optimization using mosek...\n');
  end

  if isfield(param,'Display');
    param = rmfield(param,'Display');
  end

  %ticid = tic;
  [r,res]=mosekopt(['minimize ' quiet],prob,param);
  %fprintf('mosekopt time: %g\n',toc(ticid));
  report_mosek_error(r,res);
  % extract solution from result
  x = res.sol.itr.xx;
  % set weights to solution in weight matrix
  Z = x(1:n);
  f = 0.5*sum((FF*Z).^2) + c'*Z;
end

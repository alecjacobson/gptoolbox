function x = merl_quadprog(varargin)
  % MERL_QUADPROG Optimize problems of the form:
  % 
  % min_x 1/2 x' Q x - x' h
  % 
  % s.t. x>=0
  %
  % According to the method described in "PARALLEL QUADRATIC PROGRAMMING FOR
  % IMAGE PROCESSING" [Brand & Chen 11]
  %
  % x = merl_quadprog(Q,h,'ParamterName',ParameterValue)
  %
  % Inputs:
  %   Q  n by n symmetric positive semi-definite quadratic coefficients matrix
  %   h  n by 1 linear coefficents vector (note sign)
  % Outputs:
  %   x  n by 1 solution vector
  %
  %

  Q = varargin{1};
  h = varargin{2};

  n = size(Q,1);
  assert(size(Q,2) == n);
  assert(numel(h) == n);
  h = h(:);

  max_iter = inf;
  tol = 1e-7;

  % initial guess
  x0 = ones(n,1);

  v = 3;
  while v <= numel(varargin)
    switch varargin{v}
    case 'X0'
      assert((v+1)<=numel(varargin));
      v = v+1;
      x0 = varargin{v};
    otherwise
      error(['Unknown parameter: ' varargin{v}]);
    end
    v = v+1;
  end

  %% "some r ∈ R^n ≥ 0"
  %r = sparse(n,1);
  % ri = max(Qii , ∑ j Q−i j )
  r = max(diag(Q),max(-Q,[],2));

  % "similarly" ... "si>0"
  s = ones(n,1);

  Qp = max(Q,0) + diag(r);
  Qm = max(-Q,0) + diag(r);
  assert(issparse(Qp));
  hp = max(h,0) + s;
  hm = max(-h,0) + s;

  x = x0;
  it = 0;
  while true
    x_prev = x;
    x = x.*(hp+Qm*x)./(hm+Qp*x);
    it = it + 1;
    diff = max(abs(x-x_prev));
    %fprintf('%d: %g\n',it,diff);
    if it >= max_iter
      return;
    end
    if diff<=tol
      return;
    end
  end


end

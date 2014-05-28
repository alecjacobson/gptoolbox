function [x,y,AFUN] = schur_complement(AFUN,B,C,f,g,varargin)
  % SCHUR_COMPLEMENT Solve a linear KTK system of the form:
  % [A B' * [x  = [f
  %  B C]    y]    g]
  % using the Schur complement method as described in 
  % "Numerical solution of saddle point problems" [Benzi et al. 2005]
  %
  % [x,y] = schur_complement(A,B,C,f,g)
  % [x,y] = schur_complement(AFUN,B,C,f,g)
  % 
  % Inputs:
  %   A  n by n square, positive definite matrix (typically quadratic
  %     coefficients)
  %     or
  %   AFUN  function AFUN(b) solving x = A \ b
  %   B  m by n  matrix of full row rank (typically constraints)
  %   C  m by m  matrix (typically zero)
  %   f  n list of upper right-hand side values
  %   g  m list of lower right-hand side values
  %   Optional:
  %     'Iterative' followed by whether to use an iterative solver rather than
  %       cholesky for the schur complement solve {false}.
  %     'Preconditioner' followed by a preconditioner {[]}
  %     'Dense' followed by whether to use a dense solver for the schur {true}
  %     complement solve.
  % Outputs:
  %   x  n vector of primary variable solution values
  %   y  m vector of contraint variable solution values
  %   AFUN  function AFUN(b) solving x = A \ b
  %

  function [AFUN,p] = chol_fun(A)
    [L,p,Q] = chol(A,'lower');
    % precompute transposes
    U = L';
    P = Q';
    AFUN = @(b) Q * (U \ (L \ ( P * b)));
  end

  function [AFUN,p] = ldl_fun(A)
    [L,D,Q] = ldl(A);
    % precompute transposes
    U = L';
    P = Q';
    p = sum(abs(diag(D))/max(abs(diag(D)))<1e-15);
    AFUN = @(b) Q * (U \ (D\(L \ ( P * b))));
  end

  if ~isa(AFUN,'function_handle')
    A = AFUN;
    [AFUN,p] = chol_fun(A);
    if p~=0
      warning('A is not positive definite');
      % A is not positive definite, so let's try ldl
      [AFUN,p] = ldl_fun(A);
      if p~=0
        warning('A is probably singular');
      end
    end
    [x,y] = schur_complement(AFUN,B,C,f,g,varargin{:});
    return;
  end

  if isempty(C)
    C = sparse(size(B,1),size(B,1));
  end

  % default values
  iterative = false;
  dense = true;
  preconditioner = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Iterative','Preconditioner','Dense'}, ...
    {'iterative','preconditioner','dense'});
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
  % A*x + B'*y = f
  % B*x + C*y = g
  %
  % A*x + B'*y = f
  %   multiply by B*A\
  % B*A\A*x + B*A\B'*y = f
  % B*x + B*A\B'*y = B*A\f
  if isempty(f)
    BAif = sparse(size(B,1),size(g,2));
    f = sparse(size(B,2),size(g,2));
  else
    BAif = B * AFUN(f);
  end
  %   substitute B*x = g - Cy
  % g - Cy + B*A\B'*y = B*A\f
  % (B*A\B' - C) *y = B*A\f - g
  %
  if iterative
    % computes B*A\B'*y-Cy
    BABTCFUN = @(y) B*(AFUN(B'*y))-C*y;
    tol = 1e-14;
    max_iter = 10000;
    % These seem equivalent in terms of speed. But `minres` can handles
    % symmetric indefinite matrices where `pcg` cannot.
    %y = minres(BABTFUN,BAif-g,tol,max_iter);
    y = pcg(BABTCFUN,BAif-g,tol,max_iter,preconditioner);
  else
    %% This is the bottleneck. Probably because A\B' is dense.
    BABTC = B*AFUN(B') - C;
    if dense
      BABTC = full(BABTC);
    end
    y = BABTC\(BAif-g);
    %BABTCFUN = chol_fun(BABTC);
    %assert(p==0,'A is not positive definite');
    %% y = (B*A\B'-C)\(B*A\f - g)
    %y = BABTCFUN(BAif - g);
  end
  %   substitute y in top
  % A*x + B'*((B*A\B')\(B*A\f - g)) = f
  % x = A\(f-B'*((B*A\B')\(B*A\f - g)))
  x = AFUN(f-B'*y);

end

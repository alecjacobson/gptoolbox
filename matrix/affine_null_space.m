function [N,x0] = affine_null_space(A,b,varargin)
  % Given a system Ax = b, determine a matrix N spanning the right null space
  % of A and a feasible solution x0 so that:
  %
  %     A * (N * y + x0) = b  for any y
  %
  % Inputs:
  %   A  #A by #A (sparse) matrix. Note: I'm pretty sure A must be symmetric
  %      positive semi-definite.
  %   b  #A by #b right-hand side
  %   Options:
  %     'Tol'  followed by tolerance for determine rank (what's considered zero?)
  %     'Method'  followed by either:
  %        'qr'  use QR decomposition
  %        'luq'  use LUQ decomposition
  % Outputs:
  %   N  #A by #N matrix spanning null space, where #N = #A - rank(A)
  %   x0  #A by #b, so that columns are feasible solutions
  % 

  tol = [];
  method = 'qr';
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method','Tol'},{'method','tol'});
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
    tol = max(max(size(A)) * norm(A,1) * eps,100*eps);
  end

  switch method
  case 'luq'
    % Special sparse LUQ decomposition
    [L,U,Q] = luq(A,1,tol);
    % Rank
    nc = find(any(U>tol,2),1,'last');
    if isempty(nc)
      nc = 0;
      m = size(A,1)-nc;
      x0 = ones(m,1);
      N = speye(size(A,1),m);
    else
      x0 = Q\[(U(1:nc,1:nc)\(speye(nc,size(L,1))*(L\b)));zeros(size(Q,1)-nc,1)];
      QQ = Q^-1;
      N = QQ(:,nc+1:end);
    end
    %big = max(abs(N));
    %if big > 0
    %  N = N/big;
    %end
  case 'qr'
    [Q,R,E] = qr(A');
    % Rank of A
    nc = find(any(R>tol,2),1,'last');
    % Q = [Q₁,N]
    Q1 = Q(:,1:nc);
    % A possibly non-unique solution
    x0 = Q1*(R(1:nc,1:nc)'\(E(:,1:nc)'*(b)));
    N = Q(:,nc+1:end);
  end
  % Zap anything below tolerance
  N(abs(N) < tol) = 0;
  %assert(max(abs(A*(N*rand(size(N,2),size(b,2)) + x0) - b)) < 1e-10, ...
  %  'Should span solutions to A x = b');
end

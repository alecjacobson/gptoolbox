function [N,x0] = affine_null_space(A,b,varargin)
  % Given a system Ax = b, determine a matrix N spanning the right null space
  % of A and a feasible solution x0 so that:
  %
  %     A * (N * y + x0) = b  for any y
  %
  % Inputs:
  %   A  m by n (sparse) matrix. 
  %   b  m by #b right-hand side
  %   Options:
  %     'Tol'  followed by tolerance for determine rank (what's considered
  %       zero?)
  %     'Method'  followed by either:
  %        {'qr'}  use QR decomposition of A' (robust, best understood, slowest)
  %        'luq'  use LUQ decomposition 
  %        'rq'  use QR decomposition of A (good when m << n)
  %        'svd'  use SVD decompostion (only use for small/dense matrices)
  %        'rrlu'  simplifed version of LUQ (**broken**)
  % Outputs:
  %   N  n by #N matrix spanning null space, where #N = m - rowrank(A)
  %   x0 n by #b, so that columns are feasible solutions
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
  if nargin<2 || isempty(b)
    b = zeros(size(A,1),1);
  end

  if isempty(tol)
    tol = max(max(size(A)) * norm(A,1) * eps,100*eps);
  end

  switch method
  case 'luq'
    % Special sparse LUQ decomposition
    [L,U,Q] = luq(A,1,tol);
    % Rank
    nc = find(any(abs(U)>tol,2),1,'last');
    if isempty(nc)
      nc = 0;
      m = size(A,1)-nc;
      if nargout>=2
        x0 = ones(m,1);
      end
      N = speye(size(A,1),m);
    else
      if nargout>=2
        y0 = U(1:nc,1:nc)\(speye(nc,size(L,1))*(L\b));
        x0 = Q\[y0;zeros(size(Q,1)-nc,size(b,2))];
      end
      QQ = Q^-1;
      N = QQ(:,nc+1:end);
    end
    %big = max(abs(N));
    %if big > 0
    %  N = N/big;
    %end
  case 'rrlu'
    % "Strong rank revealing LU factorizations" [Miranian & Gu 2002]
    % https://math.berkeley.edu/~luiza/RRLU.pdf
    %
    % Seems there is a typo. They write: "Nr = [-A11\A12;I_{m-k,n-k}]" But so
    % that A*Nr makes sense, I think it should be "Nr = [-A11\A12;I_{n-k,n-k}]"
    %
    %
    m = size(A,1);
    n = size(A,2);
    [L,U,P,Q] = lu(sparse(A),tol);

    NZ = find((any(abs(U)>tol,2)));
    % LUQ just checks the diagonal:
    NZluq = find(abs(diag(U))>tol);
    if ~isempty(setxor(NZ,NZluq))
      size(NZ)
      size(NZluq)
      warning('Not handling non-zero bottom right corner (use luq)');
    end

    Z = find(~(any(abs(U)>tol,2)));
    R = sparse((1:size(U,1))',[NZ;Z],1);

    C = blkdiag(R',speye(size(U,2)-size(R,2)));
    % L U               = P A Q
    % L R' R U C C'     = P A Q
    % L R' R U C        = P A Q C 
    % L R' [U11 U12; 0] = P A Q C 
    %
    % let M = [-U11\U12;I]
    %
    % [U11 U12] M = 0
    %
    % L R' [U11 U12] M = P A Q C M
    % L R' 0           = P A Q C M
    %
    % Implies
    %
    % N = Q C M --> A N = 0
    %
    UU = R*U*C;
    nc = sum((any(abs(UU)>tol,2)));
    U11 = UU(1:nc,1:nc);
    U12 = UU(1:nc,nc+1:end);
    M = [-U11\U12;speye(size(U12,2),size(U12,2))];
    N = Q*C*M;

    % We have:
    %
    % L R' [U11 U12; 0] = P A Q C 
    %
    % A x = b
    % Let x = Q C y
    % P A Q C y = P b
    % L R' [U11 U12; 0] y = P b
    % Assume: (L R') is invertible:
    % [U11 U12; 0] y = (L R')\(P b)
    % [U11 U12; 0] y = [b1;b2];
    % Invertible iff b2 = 0
    % [U11 U12; 0] y = [b1;0];
    % y = [U11\b1;y2] for any y2
    % And might as well set y2 to 0
    % y = [U11\b1; 0]
    % 
    LL = L*R';
    bb = LL\(P*b);
    y = [U11\(bb(1:nc,:));zeros(size(U,2)-nc,size(bb,2))];
    x0 = Q*C*y;

  case 'qr'
    [Q,R,E] = qr(A');
    % Rank of A
    nc = find(any(abs(R)>tol,2),1,'last');
    % Q = [Q₁,N]
    Q1 = Q(:,1:nc);
    % A possibly non-unique solution
    if nargout>=2
      x0 = Q1*(R(1:nc,1:nc)'\(E(:,1:nc)'*(b)));
    end
    N = Q(:,nc+1:end);
  case 'rq'
    % http://mathoverflow.net/a/253997/23064
    [Q,R,E] = qr(A);
    nc = find(any(abs(R)>tol,2),1,'last');
    R1 = R(1:nc,1:nc);
    R2 = R(1:nc,nc+1:end);
    n = size(A,2);
    N = E*[-(R1\R2);speye(n-nc,n-nc)];
    %assert(nargout <= 1 && 'x0 not supported for rq');
    b1 = Q(:,1:nc)'*b;
    x0 = E*[R1\b1;zeros(size(E,2)-size(R1,1),size(b1,2))];
  case 'svd'
    [U,S,V] = svd(full(A));
    % Carefully extract diagonal of S
    Sdiag = S(sub2ind(size(S),1:min(size(S)),1:min(size(S))));
    Z = abs(Sdiag)<tol;
    N = V(:,setdiff(1:end,find(~Z)));
    Sdiag(Z) = 0;
    Sdiag(~Z) = 1./Sdiag(~Z);
    % Place back into S without changing size of S
    S(sub2ind(size(S),1:min(size(S)),1:min(size(S)))) = Sdiag;
    Apinv = V*S'*U';
    x0 = Apinv * b;
  end
  % Zap anything below tolerance
  %N(abs(N) < tol) = 0
  [NI,NJ,NV] = find(N);
  N = sparse(NI,NJ,(NV>tol).*NV,size(N,1),size(N,2));
  
  %assert(max(abs(A*(N*rand(size(N,2),size(b,2)) + x0) - b)) < 1e-10, ...
  %  'Should span solutions to A x = b');
  if nargout>1 && ~(all(all(abs(A*x0-b)<tol)))
    % Should check that constraint right-hand sides are compatible:
    %   [Q,R,E] = qr(A'); 
    %   rank_A = find(any(abs(R)>tol,2),1,'last');
    %   [Q,R,E] = qr([A b]'); 
    %   rank_Ab = find(any(abs(R)>tol,2),1,'last');
    %   assert(rank_Ab <= rank_A);
    % 
    warning('MATLAB:singularMatrix', ...
      'A*x0 ~= b, this may mean that  A*x = b is impossible to satisfy');
  end
end

function [Z,E,ZZ] = lexmin(H,ff,varargin)
  % LEXMIN Solve the multi-objective minimization problem:
  %
  %     min  {E1(x), E2(x), ... , Ek(x)}
  %      x
  %
  % where
  %
  %     Ei = ½ x' H{i} x + x' f{i}
  % 
  % and Ei is deemed "more important" than Ei+1 (lexicographical ordering): 
  % https://en.wikipedia.org/wiki/Multi-objective_optimization#A_priori_methods
  %
  % Z = lexmin(H,f)
  % Z = lexmin(H,f,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   H  k-cell list of n by n sparse matrices, so that H{i} contains the
  %     quadratic coefficients of the ith energy.
  %   f  k-cell list of n by 1 vectors, so that f{i} contains the linear
  %     coefficients of the ith energy
  %     Optional:
  %       'Eps1'  followed by Tikinov regularization epsilon for each energy
  %         {0}, 1e-18 worked well once for mosekopt
  %       'Eps2'  followed by Tikinov regularization epsilon for constraints
  %         {0}, 1e-10 worked well once for mosekopt
  %       'Method' followed by one of the following method names:
  %          'fmincon'
  %          'lagrange-multiplier'
  %          {'null-space'}
  %       'NullSpaceMethod'  when using 'null-space', followed by either
  %          'qr'  use qr decomposition: better for many, big null spaces
  %          'luq'  use "luq" decomposition: better for few, small null spaces
  %       'fminconOptions' followed by cell of options to fminconv (only
  %         applicable when using 'Method','fmincon')
  %       'Z0' followed by an n by 1 initial guess
  % Outputs:
  %   Z  n by 1 solution vector
  %   E  k by 1 objective (outcome) vector
  %

  function fi = f(i)
    if iscell(ff)
      fi = ff{i};
    else
      fi = permute(ff(i,:,:),[2 3 1]);
    end
  end

  E = [];
  % Function handles used by fmincon
  function [fval,grad] = objective(X,i)
    fval = 0.5*X'*H{i}*X + X'*f(i);
    grad = H{i}*X + f(i);
  end
  function [ineqval,eqval,ineqgrad,eqgrad] = constraints(X,i)
    if i == 1
      ineqval = [];
      ineqgrad = [];
      eqval = [];
      eqgrad = [];
      return;
    end
    [ineqval,eqval,ineqgrad,eqgrad] = constraints(X,i-1);
    [fval,grad] = objective(X,i-1);
    % values in the rows
    ineqval  = [ineqval;fval-E(i-1)];
    % grads in the columns
    ineqgrad = [ineqgrad grad];
  end
  function hess = hessian(X,lambda,i)
    hess = H{i};
    for j = 2:i
      hess = hess + lambda.ineqnonlin(j-1)*H{j-1};
    end
  end

  Z = [];
  normalize = false;
  method = 'null-space';
  null_space_method = 'luq';
  fmincon_options = {};
  eps1 = 0;
  eps2 = 0;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Eps1','Eps2','Method','Normalize','NullSpaceMethod','fminconOptions','Z0'}, ...
    {'eps1','eps2','method','normalize','null_space_method','fmincon_options','Z'});
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

  % Number of energies
  k = numel(H);

  if normalize
    for i = 1:k
      si = max(abs(H{i}(:)));
      H{i} = H{i}/si;
      if iscell(ff)
        ff{i} = ff{i}/si;
      else
        ff(i,:,:) = ff(i,:,:)/si;
      end
    end
  end


  % Dimension of search space
  n = size(H{1},1);
  switch method
  case 'mosekopt'
    for i = 1:k
      clear prob;
      prob.c = f(i);
      H{i} = H{i};
      %force_psd = true;
      %if force_psd
      %  [V,S] = eig(full(H{i}));
      %  U = V;
      %  S = diag(S);
      %  S(S<0) = 0;
      %  S = diag(S);
      %  H{i} = U*S*V';
      %end
      [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(H{i} + eps1*speye(size(H{i}))));
      prob.a = sparse(0,n);
      prob.buc = [];
      prob.qcsubk = [];
      prob.qcsubi = [];
      prob.qcsubj = [];
      prob.qcval = [];
      prob.display = 'off';
      for j = 1:i-1
        [I,J,V] = find(tril(H{j}+eps2*speye(size(H{i},2))));
        prob.qcsubk = [prob.qcsubk;j*ones(numel(I),1)];
        prob.qcsubi = [prob.qcsubi;I];
        prob.qcsubj = [prob.qcsubj;J];
        prob.qcval = [prob.qcval;V];
        prob.a = [prob.a;f(j)'];
        prob.buc = [prob.buc;E(j)+1e-6];
      end
      param = [];
      [r,res] = mosekopt('minimize echo(0)',prob,param);
      if ~isempty(res.rmsg)
        error(res.rmsg);
      end
      E(i) = res.sol.itr.pobjval;
      Z = res.sol.itr.xx;
    end
  case 'fmincon'
    if isempty(Z)
      Z = -f(1);
    end
    Q = cell(k,1);
    C = cell(k+1,1);
    C{1} = @(X) [];
    GC{1} = @(X) [];
    HC{1} = @(l) 0;
    for i = 1:k
      %fprintf('i: %d\n',i);
      Q{i} = @(X) 0.5*(X'*(H{i}*X)) + X'*f(i);
      G{i} = @(X) (H{i}*X) + f(i);
      options = optimoptions( ...
        'fmincon', ...
        'Algorithm','interior-point', ...
        'Display', 'off', ...
        'GradObj','on','Hessian','on', ...
        'GradConstr','on', ...
        'Hessfcn',@(X,lambda) hessian(X,lambda,i), ...
        fmincon_options{:});
      % It's important that quadratic constraints from the previous
      % optimizations are specified as *in*-equality constraints. They're
      % approximate and conservative, so forcing equality is
      % counter-productive.
      warning('off');
      [Z,E(i)] = fmincon( ...
         @(X) objective(X,i), ...
         Z, ...
         [],[], ... % A,B
         [],[], ... % Aeq,Beq
         [],[], ... % LB,UB
         @(X) constraints(X,i), ...
         options);
      warning('on');
      % no early exit
    end
  case 'lagrange-multiplier'
    % This method is bad for multiple reasons: 
    %   1. qr factorization on the ith set of constraints will be very
    %      expensive,
    %   2. this is due to the fact that the number of non-zeros will be
    %      O(n2^i), and
    %   3. the eventual solve will not necessarily behave well because the
    %      constraints are not full rank (!duh!).
    C = H{1};
    D = -f(1);
    for i = 1:k
      [Q,R,E] = qr(C');
      nc = find(any(R,2),1,'last');
      if nc == size(C,1) || i==k
        Z = C \ D;
        Z = Z(1:n);
        break;
      end
      %     min               ½x'Hix + x'fi + 0*λ₁ + 0*λ₂ + ... + 0*λi-1
      %      x,λ₁,λ₂,...λi-1
      %      s.t.             Ci-1 * [x; λ₁; λ₂; ...; λi=1] = Di-1
      % or
      %     min               ½[x;λ₁;...;λi-1]'A[x;λ₁;...;λi-1] - [x;λ₁;...;λi-1]'B
      %      x,λ₁,λ₂,...λi-1                                    ^
      %                                                         |------------negative
      %     s.t.              Ci-1 * [x; λ₁; λ₂; ...; λi-1] = Di-1
      A = sparse(size(C,1),size(C,1));
      A(1:n,1:n) = H{i+1};
      B = zeros(size(C,1),1);
      B(1:n) = -f(i+1);
      C = [A C';C sparse(size(C,2),size(C,2))];
      D = [B;D];
    end
  case 'null-space'
    % Start with "full" search space and 0s as feasible solution
    N = 1;% aka speye(n,n)
    Z = zeros(n,1);
    for i = 1:k
      % Original ith energy: ½ x' Hi x + x' fi
      Hi = H{i};
      fi = f(i);
      % Restrict to running affine subspace, by abuse of notation:
      %     x = N*y + z
      fi = N'*(Hi*Z+fi);
      Hi = N'*Hi*N;
      [Ni,Y] = affine_null_space(Hi+eps2*speye(size(Hi)),-fi,'Method',null_space_method);
      % Update feasible solution
      Z = N*Y + Z;
      if size(Ni,2) == 0
        % Z is full determined, exit loop early
        break;
      end
      % Otherwise, N spans the null space of Hi
      N = N*Ni;

      %Z = N*Y + Z;
      %if size(Ni,2) == 0
      %  break;
      %end
      %N = N*Ni;

    end
    % (If i<k then) the feasible solution Z is now the unique solution.
  case 'kanoun'
    % As described in [De Lasa et al. 2010] Section 8.1
    ZZ = cell(k,1);
    Aeq = [];
    Beq = [];
    for i = 1:k
      % I'll use j where De Lasa uses k
      %ZZ{i} = min_quad_with_fixed(0.5*H{i},f(i),[],[],Aeq,Beq);
      ZZ{i} = quadprog(H{i},f(i),[],[],Aeq,Beq);
      if isempty(ZZ{i}) || size(ZZ{i},1) < size(H{i},1)
        Z = ZZ{i-1};
        break;
      end
      % H{i} = - S*P*L*D*L'*P'*S'
      % √H{i} = sqrt(D)*L'*P'*S' 
      [L,D,P,S] = ldl(H{i}+eps2*speye(size(H{i})));
      D = diag(D);
      D(D<0) = 0;
      nc = find(D,1,'last');
      D = diag(D);
      % only keep those with non-zero D
      D = D(1:nc,:);
      Ai = sqrt(D)*L'*P'*S';
      Aeq = [Aeq;Ai];
      Beq = [Beq;Ai*ZZ{i}];
      E = zeros(i,1);
      for j = 1:i
        E(j) = 0.5*(ZZ{i}'*(H{j}*ZZ{i})) + ZZ{i}'*f(j);
      end
    end
    Z = ZZ{end};
  end

  if nargout>1
    E = zeros(k,1);
    for i = 1:k
      E(i) = 0.5*(Z'*(H{i}*Z)) + Z'*f(i);
    end
  end
end

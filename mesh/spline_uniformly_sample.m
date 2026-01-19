function [U,I,t,E] = spline_uniformly_sample(P,C,n,varargin)
  % [U,I,t,E] = spline_uniformly_sample(P,C,n,varargin)
  %
  % Inputs:
  % P  #P by dim list of vertex positions
  % C  #C by 4 list of cubic bezier control points indices into P
  % n  number of samples
  %   Optional:
  %     'Fsum'  followed by the sum of the arc lengths of the curve
  %     'Stagger'  followed by whether to stagger the samples
  %     'Quadrature'  followed by #Q by 2 list of quadrature points and weights
  %     'NumQuadrature'  followed by number of quadrature points to use
  %     'MaxIter'  followed by maximum number of iterations
  %     'Tol'  followed by tolerance for convergence
  %     'GradTol'  followed by tolerance for convergence
  % Outputs:
  %   U  n by dim list of uniformly sampled points
  %   I  n list of indices into C
  %   t  n list of parameters into cubic bezier curves
  %   E  #E by 2 list of edges connecting neighboring samples

  stagger = [];
  f_tol = 1e-16;
  grad_tol = 1e-16;
  nq = [];
  Q = [];
  max_iter = 100000;
  given_Fsum = [];
  viz = false;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Fsum','Stagger','Quadrature','NumQuadrature','MaxIter','Tol','GradTol'}, ...
    {'given_Fsum','stagger','Q','nq','max_iter','f_tol','grad_tol'});
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


  if (numel(nq) == 1 || isempty(nq)) && isempty(Q) 
    if isempty(nq)
      nq = 10;
    end
    Q = nan(nq,2);
    [Q(:,1),Q(:,2)] = gauss_legendre_quadrature(nq);
  end
  nq = size(Q,1);
  assert(size(Q,2) == 2);

  if ~all(C(1:end-1,4) == C(2:end,1))
    % first call on each component individually
    [~,CC] = connected_components(C);
    nk = max(CC);
    if nk > 1

      S = sparse(1:size(C,1),CC,1,size(C,1),max(CC));
      l = zeros(nk,1);
      for c = 1:nk
        Ic = find(S(:,c));
        Cc = C(Ic,:);
        [Pc,~,~,Cc] = remove_unreferenced(P,Cc,true);
        l(c) = sum(spline_arc_lengths(Pc,Cc));
      end

      nn = max(ceil(n*l/sum(l)),2);
      % sum(nn) not necessarily equal to n due to rounding
      U = [];
      I = [];
      t = [];
      E = [];
      for c = 1:nk
        Jc = find(S(:,c));
        Cc = C(Jc,:);
        [Pc,~,~,Cc] = remove_unreferenced(P,Cc,true);
        [Uc,Ic,tc,Ec] = spline_uniformly_sample(Pc,Cc,nn(c),varargin{:});
        E = [E size(U,2)+Ec'];
        U = [U Uc'];
        I = [I;Jc(Ic)];
        t = 'not supported';
      end
      U = U';
      E = E';
      return;
    else
      error('Single component must be a chain in order');
    end
  end

  if isempty(stagger)
    stagger = C(1,1) == C(end,end);
  end

  on_loop = C(1,1) == C(end,end);

  if on_loop
    assert(stagger,'Loops should use ''Stagger'',true');
    C(end,end) = size(P,1)+1;
    P = [P;P(1,:)];
    [U,I,t,E] = spline_uniformly_sample(P,C,n+1,'Stagger',false);
    U = U(1:end-1,:);
    I = I(1:end-1);
    t = t(1:end-1);
    E = [E(1:end-1,:);E(end,1) E(1,1)];
    return;
  end

  % Naive initial guess
  nc = size(C,1);
  gt = initial_guess();
  if stagger && ~on_loop
    %gt = nc*linspace(1/(2*n),1-1/(2*n),n)';
    w = [0.5;ones(n-1,1);0.5];
    w = w./sum(w);
  else
    %gt = linspace(0,1*nc,n)';
  end
  g2I = @(gt) min(floor(gt)+1,nc);
  g2l = @(gt,I) deal(I,gt-(I-1));
  g2l = @(gt) g2l(gt,g2I(gt));
  l2g = @(I,l) I-1+l;


  for iter = 1:max_iter
    if viz
      U = [];
      [I,lt] = g2l(gt);
      cumlens = cumsum([0;accumarray(I,1,[nc 1])]);
      for c = 1:nc
        if cumlens(c+1)-cumlens(c) == 0
          continue;
        end
        ti = lt(cumlens(c)+1:cumlens(c+1));
        U = [U;cubic_eval(P(C(c,:),:),ti)];
      end
      [vV,vE] = spline_to_poly(P,C,0.01);
      tsurf(vE,vV);
      hold on;
      sct(P(C(:,[1 4]),:),'r','filled');
      sct(U,'b','filled');
      hold off;
      set(gca,'YDir','reverse');
      axis equal;
      drawnow;
      pause
    end

    [res,J] = residuals(P,C,gt);
    %f_J = fd_jacobian(@(gt) residuals(P,C,gt),gt);
    %full(J)
    %full(f_J)
    %keyboard
    %error

    dfdt = J'*res;
    if stagger
      %assert(false)
    else
      dfdt([1 end]) = 0;
    end
    method = 'gn';
    switch method
    case 'gd'
      dt = -dfdt;
    case 'gn'
      %dt = -(pinv(J,1e-10)*dfdt);
      H = J'*J;
      dt = zeros(size(dfdt));
      if stagger
        dt = -(H\dfdt);
      else
        dt(2:end-1) = -(H(2:end-1,2:end-1)\dfdt(2:end-1));
      end
    end

    f = @(gt) sum(residuals(P,C,gt).^2,'all');
    [ss,gt,f_gt] = backtracking_line_search(f,gt,dfdt,dt,0.01,0.5);
    gt = sort(gt);
    if f_gt<f_tol
      break;
    end
    if ss == 0
      error('Line search failed');
    end
      %fprintf('%d: %g %g %g %g\n',iter,ss,f_gt,norm(dfdt,inf),norm(dt,inf));


  end

  U = [];
  [I,lt] = g2l(gt);
  cumlens = cumsum([0;accumarray(I,1,[nc 1])]);
  for c = 1:nc
    ti = lt(cumlens(c)+1:cumlens(c+1));
    U = [U;cubic_eval(P(C(c,:),:),ti)];
  end
  t = lt;
  if on_loop
    E = [1:size(U,1);2:size(U,1) 1]';
  else
    E = [1:size(U,1)-1;2:size(U,1)]';
  end

  function [res,J] = residuals(P,C,gt)
    if on_loop
      % wrap around
      gt(gt>nc) = gt(gt>nc)-nc;
      gt(gt<0) = gt(gt<0)+nc;
    end
    gt = sort(gt);

    compute_J = nargout > 1;
    if ~compute_J
      if any(gt>size(C,1)) || any(gt<0)
        res = inf;
        return;
      end
    end

    [I,lt] = g2l(gt);
    cumlens = cumsum([0;accumarray(I,1,[nc 1])]);

    if stagger
      if on_loop
        F = zeros(n,1);
      else
        F = zeros(n+1,1);
      end
    else
      F =   zeros(n-1,1);
    end


    if compute_J
      JI = [];
      JJ = [];
      JV = [];
    end
    Fsum = 0;
    for c = 1:nc
      ti = lt(cumlens(c)+1:cumlens(c+1));
      si = [0;ti;1];
      if compute_J
        [Fi,dFida,dFidb] = cubic_arc_length(P(C(c,:),:),Q,si(1:end-1),si(2:end));
      else
        [Fi] = cubic_arc_length(P(C(c,:),:),Q,si(1:end-1),si(2:end));
      end
      Fsum = Fsum + sum(Fi);
      into = 0:numel(Fi)-1;
      keep = 1:numel(Fi);
      if on_loop
        % There's a smart way to fix this for loops but I'm confused. Loops are
        % handled above by faking them as non-looping.
        assert(false)
      end

      if stagger
        from = into;
        into = into+1;
      else
        if c == 1
          into = into(2:end);
          keep = keep(2:end);
        end
        if c == nc
          into = into(1:end-1);
          keep = keep(1:end-1);
        end
        from = into;
      end

      FI = cumlens(c)+into';
      tJ = cumlens(c)+from(2:end)';

      if compute_J
        JIa = FI(2:end);
        JJa = tJ;
        JVa = dFida(keep(2:end));
        JIb = FI(1:end-1);
        JJb = tJ;
        JVb = dFidb(keep(1:end-1));
        JI = [JI(:);JIa(:);JIb(:)];
        JJ = [JJ(:);JJa(:);JJb(:)];
        JV = [JV(:);JVa(:);JVb(:)];
      end

      %F(FI) = F(cumlens(c)+into') + Fi(keep);
      F(FI) = F(FI) + Fi(keep);

      if compute_J
        J = sparse(JI,JJ,JV,size(F,1),n);
      end
    end
    if ~isempty(given_Fsum)
      Fsum = given_Fsum;
    end
    if stagger && ~on_loop
      h = Fsum.*w;
    else
      h = Fsum/numel(F);
    end
    res = F-h;
  end

  function gt = initial_guess()
    l = spline_arc_lengths(P,C);
    total_l = sum(l);
    p = l/total_l;
    X = round(n.*p);
    %[~,I] = sort( abs(X/n -p) );
    while true
      deficit = n - sum(X);
      if deficit == 0
        break;
      end
      % when deficit > 0 then we need to add elements to X
      if deficit > 0
        [~,I] = sort([n.*p-X],'descend');
      else
        [~,I] = sort([n.*p-X],'ascend');
      end
      X(I(1:min(abs(deficit),end))) = X(I(1:min(abs(deficit),end))) + sign(deficit);
    end
    gt = [];
    for c = 1:nc
      xc = X(c);
      lt = linspace(1/(2*xc),1-1/(2*xc),xc)';
      gt = [gt;lt + (c-1)];
    end
    if ~stagger && ~on_loop
      gt(1) = 0;
      gt(end) = nc;
    end
  end

end

function [t,F,U,dtdC,dFdC,dUdC] = cubic_uniformly_sample(C,n,varargin)
  % [t,F,U] = cubic_uniformly_sample(C,n)
  %
  % Inputs:
  %  C  #C by 4 list of cubic control points
  %  n  number of samples
  %    Optional:
  %      'NumQuadrature' followed by number of quadrature points {10}
  %      'MaxIter' followed by maximum number of iterations {10}
  %      'Tol' followed by tolerance for function value {1e-16}
  %      'GradTol' followed by tolerance for gradient {1e-16}
  %      'Stagger' followed by whether to compute uniform samples offset half a
  %        step from each endpoint.
  %      'Fsum'  followed by a previously computed sum of output lengths. This
  %        is a debugging option so that we can pass the sum of lengths without
  %        differentiating through it (only relevant for 'Stagger',true).
  %        Otherwise, the derivatives will include terms corresponding to
  %        changing the overall length of the curve.
  % Outputs:
  %  t  n list of parameter values
  %  F  n-1 list of segment lengths
  %  U  n by size(C,2) list of uniformly sampled points
  %
  % Example:
  %  % padded samples so that each *sample* corresponds to *center* of an equal
  %  % length segment.
  %  % (perhaps consider just using [½ 1 1 … 1 1 ½] weights instead)
  %  [t,F,P] = cubic_uniformly_sample(C,2*n+1);
  %  t = t(2:2:end);
  %  F = F(1:2:end) + F(2:2:end);
  %  P = P(2:2:end,:);

  % number of quadrature points
  f_tol = [];
  grad_tol = 1e-16;
  nq = [];
  Q = [];
  max_iter = 20;
  stagger = false;
  Fsum = [];
  t = [];


  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Fsum','Stagger','Quadrature','NumQuadrature','MaxIter','Tol','GradTol','InitialGuess'}, ...
    {'Fsum','stagger','Q','nq','max_iter','f_tol','grad_tol','t'});
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
  if stagger
    assert(n>=1);
  else
    assert(n>=2);
  end

  if isempty(f_tol)
    bbd = norm(max(C)-min(C));
    f_tol = 1e-16*bbd;
  end

  if isempty(t)
    if stagger
      t = linspace(1/(2*n),1-1/(2*n),n)';
    else
      t = linspace(0,1,n)';
    end
  end

  if stagger
    s_fun = @(t) [0;t;1];
    ns = n+2;
  else
    s_fun = @(t) t;
    ns = n;
  end
  JI = repmat((1:ns-1),2,1)';
  JJ = [1:ns-1;2:ns]';
  w = [0.5;ones(n-1,1);0.5];
  w = w./sum(w);

  %if isempty(Fsum)
  %  Fsum = cubic_arc_length(C,nq*(ns-1),0,1);
  %end


  if (numel(nq) == 1 || isempty(nq)) && isempty(Q) 
    if isempty(nq)
      nq = 10;
    end
    Q = nan(nq,2);
    [Q(:,1),Q(:,2)] = gauss_legendre_quadrature(nq);
  end
  nq = size(Q,1);
  assert(size(Q,2) == 2);

  for iter = 1:max_iter
    s = s_fun(t);
    [F,dFda,dFdb] = cubic_arc_length(C,Q,s(1:end-1),s(2:end));
    assert(numel(F) == ns-1);

    J = sparse(JI,JJ,[dFda dFdb],ns-1,ns);

    % Making Fsum dependent on C is of course captured by numerical
    % differentiation and may result in unexpected differences from assuming
    % that Fsum is constant. (especially in the 'Stagger',true case). When
    % comparing to finite differences for example, it's better to precompute
    % Fsum outside of the differentiated call.
    %
    % It doesn't seem necessary to recompute this sum every iteration.
    if isempty(Fsum)
      Fsum = sum(F);
    end 
    if stagger
      h = Fsum.*w;
    else
      h = Fsum/(ns-1);
    end
    dfdt = J'*(F-h);
    if stagger
      dfdt = dfdt(2:end-1);
    end

    f = @(s) 0.5*sum((cubic_arc_length(C,Q,s(1:end-1),s(2:end)) - h).^2,'all');
    f = @(t) f(s_fun(t));

    if norm(dfdt,inf) < grad_tol
      %warning('converged (iter=%d)',iter);
      break;
    end

    method = 'gd';
    method = 'sobolev';
    method = 'gn';
    switch method
    case 'gd'
      dt = -dfdt;
    case 'sobolev'
      assert(~stagger);
      dt([1 end]) = 0;
      E = [1:numel(t)-1;2:numel(t)]';
      A = adjacency_matrix(E);
      L = diag(sum(A,2))-A;
      dt = [0;-(L(2:end-1,2:end-1) \ dfdt(2:end-1));0];
    case 'gn'
      H = J'*J;
      if stagger
        dt = -(H(2:end-1,2:end-1) \ dfdt);
      else
        dt = [0;-(H(2:end-1,2:end-1) \ dfdt(2:end-1));0];
      end
    end

    [ss,t,ft] = backtracking_line_search(f,t,dfdt,dt,0.3,0.5);
    if nargout<=3
      %fprintf('  %d: %g %g %g %g\n',iter,ss,ft,norm(dfdt,inf),norm(dt,inf));
    end

    %clf;
    %plot_cubic(C,[],[],{{'LineWidth',1},{'LineWidth',0.5}});
    %hold on;
    %sct(cubic_eval(C,t),'or','filled');
    %hold off;
    %axis equal;
    %fprintf('iter %d: %g %g\n',iter,s,f(t));
    %drawnow;

    if ft<f_tol
      % recompute after step one last time
      s = s_fun(t);
      F = cubic_arc_length(C,Q,s(1:end-1),s(2:end));
      break;
    end
    if ss == 0
      error('line search failed');
    end
  end

  if nargout>2
    U = cubic_eval(C,t);

    if nargout>3
      s = s_fun(t);
      [F,dFda,dFdb,dFdC] = cubic_arc_length(C,Q,s(1:end-1),s(2:end));
      if stagger
        h = sum(F).*w;
      else
        h = mean(F);
      end
      dFdt = sparse(JI,JJ,[dFda dFdb],ns-1,ns);
      if stagger
        dFdt = dFdt(:,2:end-1);
      end
      %F_func = @(s) cubic_arc_length(C,Q,s(1:end-1),s(2:end));
      %F_func = @(t) F_func(s_fun(t));
      %f_dFdt = fd_jacobian(F_func,t);

      %full(dFdt)
      %full(f_dFdt)
      %error

      % H is |t| by |t|
      H = dFdt'*dFdt;
      % dFdC is |F| by |C|
      %f_dFdC = fd_jacobian(@(z) cubic_arc_length(reshape(z,size(C)),Q,s(1:end-1),s(2:end)),C(:));
      %max(abs(dFdC - reshape(f_dFdC,size(dFdC))),[],'all')
      % dFdSᵀ f_dFdC is |t| by |C|
      
      %% I hope I'm not getting lucky on the example I'm testing but it seems
      %% that contraction is only non-zero at the endpoints, which we're zeroing
      %% out below anyway. So just skip this (which thankfully means skipping the
      %% mixed second partials).
      %contraction = zeros(numel(t),numel(C));
      %for i = 1:numel(C)
      %  % d²F/dtdCᵢ is |F| by |t|
      %  f_ddFdtdCi = reshape( ...
      %    cs_jacobian(@(z) reshape(cubic_arc_length_jacobian(set(C,i,z),Q,t),[],1),C(i)),size(dFdt));
      %  % d²F/dtdC is |F| by |t| by |C|
      %  % d²F/dtdCᵀ ⋅ F  would be |t| by |C|
      %  % d²F/dtdCᵢᵀ ⋅ F  is |t| by 1
      %  contraction(:,i) = f_ddFdtdCi' * (F-h);
      %  %                                 ^--- this is zero...
      %end
      %rhs = (contraction + dFdt'*reshape(dFdC,[],numel(C)));
      rhs = (dFdt'*reshape(dFdC,[],numel(C)));

      if stagger
        %max(abs(contraction),[],'all')
        dtdC = -(H \ rhs);
      else
        %max(abs(contraction(2:end-1,:)),[],'all')
        rhs([1 end],:) = 0;
        dtdC = zeros(n,numel(C));
        % This involves 8 back-substitutions. It's likely that this dtdC 
        % is getting multiplied by something on its left with only one row.
        % Returning a function handle to the action of this multiplication would
        % be more efficient. Adjoint method.
        %
        % H(2:end-1,2:end-1) is nq-2 by nq-2
        % rhs(2:end-1,:) is nq-2 by 4×2=8
        dtdC(2:end-1,:) = -(H(2:end-1,2:end-1) \ rhs(2:end-1,:));
      end
      dtdC = reshape(dtdC,[n size(C)]);


      %f_dtdC = reshape( ...
      %  fd_jacobian( ...
      %    @(z) reshape( ...
      %      extract_arg(1,@cubic_uniformly_sample,reshape(z,size(C)),n,'Stagger',stagger),[],1),C(:)),[n size(C)]);
      %    d = round((rhs - -H*reshape(f_dtdC,n,[]))*10000)/10000;
      %    [n d(1) sum(F)*[t(1) 1-t(end)]]
      %    %round((rhs - -H*reshape(f_dtdC,n,[]))*10000)/10000
      %    %rhs ./ (-H*reshape(f_dtdC,n,[]))
      %    error

      % U =  ...
      %   1*(1-t).^3.*t.^0.*C(1,:) +  ...
      %   3*(1-t).^2.*t.^1.*C(2,:) +  ...
      %   3*(1-t).^1.*t.^2.*C(3,:) +  ...
      %   1*(1-t).^0.*t.^3.*C(4,:);
      b1 = 1*(1-t).^3.*t.^0;
      b2 = 3*(1-t).^2.*t.^1;
      b3 = 3*(1-t).^1.*t.^2;
      b4 = 1*(1-t).^0.*t.^3;

      dUdt = cubic_tangent(C,t);
      dUdC = zeros(n,size(C,2),size(C,1),size(C,2));
      dUdC(:,1,:,1) = cat(3,b1,b2,b3,b4);
      dUdC(:,2,:,2) = dUdC(:,1,:,1);
      % not clear how this effects adjoint method above. I guess it means dUdC
      % should be returned as a multplication action function handle.
      dUdC = dUdC + dUdt.* reshape(dtdC,[n 1 size(C)]);

      % total dFdC = ∂F/∂C + ∂F/∂t ∂t/∂C
      % dFdC is |F| by |C| by 2
      % dFdt is |F| by |t|
      % dtdC is |t| by |C| by 2
      % dFdC = dFdC + dFdt : dtdC
      dFdC = dFdC + reshape(dFdt*reshape(dtdC,numel(t),[]),size(dFdC));



    end
  end

  function C = set(C,i,v)
    C(i) = v;
  end

  function J = cubic_arc_length_jacobian(C,Q,t)
    l_s = s_fun(t);
    l_ns = numel(l_s);
    l_n = numel(t);
    l_JI = repmat((1:l_ns-1),2,1)';
    l_JJ = [1:l_ns-1;2:l_ns]';
    [~,l_dFda,l_dFdb] = cubic_arc_length(C,Q,l_s(1:end-1),l_s(2:end));
    J = sparse(l_JI,l_JJ,[l_dFda l_dFdb],l_ns-1,l_ns);
    if stagger
      J = J(:,2:end-1);
    end
  end

end

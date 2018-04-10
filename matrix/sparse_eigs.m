function [Phi,lambda] = sparse_eigs(L,D,k,mu,varargin)
  % Sparse (aka "compressed") "eigenmodes" minimizing:
  %
  % min_φ  tr( φ' L φ ) + μ ‖φ‖₁  subject to φ' D φ = I
  %
  % [Phi,lambda] = sparse_eigs(L,D,k,mu,varargin)
  %
  % Inputs:
  %   L  #L by #L sparse, semi-positive definite "stiffness" matrix 
  %   D  #L by #L diagonal positive definite "mass" matrix
  %   k  number of modes requestedd
  %   mu  sparsity parameter (larger is more sparse) {10}
  %   Optional:
  %     'Method'  followed by one of the following:
  %       'neumann'  This ADMM optimization is implemented according to
  %          "Compressed Manifold Modes for Mesh Processing" [Neumann et al.
  %          2014]. This method only works for {'Constrained',false}
  %       'brandt'  This follows the implementation in "Compressed Vibration
  %         Modes of Elastic Bodies" [Brandt & Hildebrandt 2017]
  %     'Constrained'  followed by whether to solve the L1 "constrained" problem
  %        where for all i>1, instead of minimizing μ‖φi‖₁, we add constrain its
  %        L1 norm to match that of the first mode (i=1): ‖φi‖₁ = ‖φ₁‖₁
  %        "Compressed Vibration Modes of Elastic Bodies" [Brandt & Hildebrandt
  %        2017]
  % Outputs:
  %   Phi  #L by k modes sorted in decreasing order by energy value
  %   lambda  k by k diagonal matrix of energy values (in some way analogous to
  %     eigen values, though L*Phi ≠ lambda*Phi)
  %

  method = 'neumann';
  constrained = false;
  max_iter = 2000;
  % Whether to use admm.m, augh... it's slower because of function handles!
  admm_function = false;
  quadprog_max_iter = 100;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'ADMM','Method' ,'Constrained','MaxIter','QuadProgMaxIter'}, ...
    {'admm_function','method','constrained','max_iter','quadprog_max_iter'});
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

  % PRECOMPUTATION
  switch method
  case 'neumann'
    n = size(L,1);
    assert(isdiag(D));
    Dsqrt = diag(sqrt(diag(D)));
    isqrt = @(D) diag(sqrt(diag(D)).^-1);
    Disqrt = isqrt(D);
    I = speye(n);
    num_refactors = 0;
    shrink = @(u,delta) sign(u).*max(abs(u)-delta,0);
    %M = speye(size(D,1));
    M = D;
  end

  % Some functions for the ADMM 'neumann' method
  function [Phi,data] = argmin_Phi(ES,UEUS,rho,data)
    E = ES(1:n,:);
    S = ES(n+(1:n),:);
    UE = UEUS(1:n,:);
    US = UEUS(n+(1:n),:);
    % Based on the Python code:
    Y = Dsqrt*(E-UE + S-US);
    % Matlab:  [U,S,V] = svd(A) --> U*S*V' = A
    % Numpy:  [U,S,VT] = svd(A) --> U*S*VT = A
    [V,W,VT] = svd(Y'*Y);
    Wisqrt = isqrt(W);
    Psi = (Y * (V * Wisqrt * VT'));
    Phi = Disqrt * Psi;
  end
  function [ES,data] = argmin_ES(Phi,UEUS,rho,data)
    UE = UEUS(1:n,:);
    US = UEUS(n+(1:n),:);
    if isempty(data)
      data.rho = inf;
    end
    data.rho_prev = data.rho;
    data.rho = rho;
    % Python code updates S before E... shouldn't matter
    %% Based on paper (should be same)
    % v = Phi + US;
    % S = sign(v).*max(abs(v) - mu/rho,0);
    % Based on python code
    S = shrink(Phi + US,(mu/rho)*full(diag(M)));
  
    % Update E (11)
    if data.rho ~= data.rho_prev
      A = data.rho*I + L + L';
      [cL,p,cS] = chol(A,'lower');
      cU = cL';
      cP = cS';
      cQ = cS;
      data.chol_solve = @(X) cQ * (cU \ (cL \ ( cP * X)));
      num_refactors = num_refactors+1;
    end
    E = data.chol_solve(data.rho*(Phi + UE));
    %E = (rho*I + L + L')\( rho*(Phi + UE) );
    ES = [E;S];

    %clf;
    %hold on;
    %off = [max(VV(:,1))-min(VV(:,1)) 0 0];
    %R = axisangle2matrix([1 0 0],-pi/2);
    %tsurf(F,(VV+0*off)*R,'CData',Phi(:,1),fphong,fsoft,'EdgeColor','none');
    %tsurf(F,(VV+1*off)*R,'CData',E(:,1),fphong,fsoft,'EdgeColor','none');
    %tsurf(F,(VV+2*off)*R,'CData',S(:,1),fphong,fsoft,'EdgeColor','none');
    %hold off;
    %axis equal;
    %caxis([-1 1]*max(max(abs(Phi(:,1)))));
    %camproj('persp');
    %set(gcf,'Color',0.9*[1 1 1]);
    %set(gca,'Visible','off');
    %l = light('Position',[-0.5 10 20],'Style','local');
    %s = ...
    %  add_shadow([],l,'Fade','local','Color',get(gcf,'Color')*0.9,'BackgroundColor',get(gcf,'Color'));
    %view(0,0);
    %colorbar
    %drawnow;
  end

  switch method 
  case 'neumann'
    % Even if I measure the L1 norm of S weighted by vertex area, I'm getting bias
    % toward higher resolution areas. If I warm start from a regular mesh, this
    % almost goes away suggesting it's a local minimum issue, but there's still
    % some assymmetric bais. My only idea to fix this so far is to also measure
    % the penalties on the auxiliary variables using weighted vertex areas...
    % (S-Phi)'*M(S-Phi) instead of (S-Phi)'*(S-Phi). This should not change the
    % actual energy landscape, but maybe knock things into a more fair local
    % minimum...
    %

    if admm_function
      [Phi,~] = eigs(L,D,k,'sm');
      % Prepare initial state
      state = struct();
      state.X = Phi;
      state.Z = [Phi;Phi];
      state.U = zeros(2*n,k);
      state.rho_prev = nan;
      state.rho = 1;
      state.argmin_X_data = [];
      state.argmin_Z_data = [];
      Phi = admm( ...
        @argmin_Phi,@argmin_ES,[I;I],-blkdiag(I,I),zeros(2*n,k),state, ...
        'MaxIter',max_iter);
    else
      % From Python code.
      tol_abs = 1e-8;
      tol_rel = 1e-6;
      % Python codes uses qr factorization of random numbers in [-1,1] range,
      % presummably to start with orthogonal vectors
      %[Psi,~] = qr(2*rand(n,k)-1,0);
      %Phi = Disqrt*Psi;
      % This is (more) deterministic. The cost is negligible in comparison to the
      % following.
      [Phi,~] = eigs(L,D,k,'sm');
      vec = @(X) X(:);
      frob_sqr = @(X) vec(X)'*vec(X);
  
      S = Phi;
      Psi = Dsqrt * Phi;
      E = S;
      U = zeros(2*n,k);
      
      rho_prev = nan;
      rho = 1;
      % paper says multiply by n, but python code divides by n.
      % This does seem to be the right units. 
      %mu = mu/n;
      
      % From Python Code
      check_interval = 10;
      
      for iter = 1:max_iter
      
        % Update Phi (10)
        UE = U(1:n,:);
        US = U(n+(1:n),:);
        %% Based on the paper:
        %Y = 0.5*(S - US + E - UE);
        %A = Dsqrt*Y;
        %[V,W] = eig(A'*A);
        %% Now  V*W*V' = (Dsqrt*Y)'*(Dsqrt*Y)
        %Phi = Disqrt * (Y*V*isqrt(W)*V');
  
        % Based on the Python code:
        Y = Dsqrt*(E-UE + S-US);
        % Matlab:  [U,S,V] = svd(A) --> U*S*V' = A
        % Numpy:  [U,S,VT] = svd(A) --> U*S*VT = A
        [V,W,VT] = svd(Y'*Y);
        Phi_prev = Phi;
        Wisqrt = isqrt(W);
        Psi = (Y * (V * Wisqrt * VT'));
        Phi = Disqrt * Psi;
  
        %e10 = @(Phi,E,S,UE,US) rho/2 * frob_sqr([Phi;M*Phi] - [E;S] + [UE;US]);
        %e10_before = e10(Phi,E,S,UE,US);
      
        % Python code updates S before E... shouldn't matter
        % Update S (12)
        Sprev = S;
        %% Based on paper (should be same)
        % v = Phi + US;
        % S = sign(v).*max(abs(v) - mu/rho,0);
        % Based on python code
        S = shrink(Phi + US,(mu/rho)*full(diag(M)));

        % TODO: it'd be nice to also handle the brandt-style constrained problem
        % in this setup. We would need to:
        %
        % min_S ρ/2‖Φ - S + US‖² subject to ‖Sj‖₁ = μ ∀j
      
        % Update E (11)
        if rho ~= rho_prev
          A = rho*I + L + L';
          [cL,p,cS] = chol(A,'lower');
          cU = cL';
          cP = cS';
          cQ = cS;
          chol_solve = @(X) cQ * (cU \ (cL \ ( cP * X)));
          num_refactors = num_refactors+1;
        end
        Eprev = E;
        E = chol_solve(rho*(Phi + UE));
        %Eprev = E;
        %E = (rho*I + L + L')\( rho*(Phi + UE) );
      
        % Update U (13)
        U = U+[Phi;Phi]-[E;S];
        % update rho (3.4.1 of "Distributed Optimization and Statistical Learning via
        % the Alternating Direction Method of Multipliers" [Boyd et al. 2011]
        bmu = 5;
        btao_inc = 2;
        btao_dec = 2;
        rho_prev = rho;
        dual_residual = sqrt(rho*sum((S(:)-Sprev(:)).^2+(E(:)-Eprev(:)).^2));
        residual = sqrt(sum((Phi(:)-E(:)).^2+(Phi(:)-S(:)).^2));
        if mod(iter,check_interval) == 0
          if residual > bmu*dual_residual
            rho = btao_inc*rho;
            % From python code: https://github.com/tneumann/cmm/blob/master/cmmlib/cmm.py
            U = U/btao_inc;
          elseif dual_residual > bmu*residual
            rho = rho/btao_dec;
            U = U*btao_dec;
          end
        end
        % From Python code / This is not _quite_ what Boyd suggests...
        eps_pri = sqrt(k*2)*tol_abs + tol_rel*max(sqrt(sum([Phi;E;S].^2,2)));
        eps_dual = sqrt(k)*tol_abs + tol_rel*rho*sqrt(sum(U(:).^2));
        if residual < eps_pri && dual_residual < eps_dual
          break;
        end
      end
    end
  case 'brandt'
    %mu = 1e2;
    n = size(L,1);
    assert(n == size(L,2));
    assert(n == size(D,1));
    assert(n == size(D,2));
    % use eigs to deterministically choose initial vectors, but avoid first one
    % (constant for the Laplacian)
    [U,~] = eigs(L,D,k+1,'sm');
    U = abs(fliplr(U(:,1:k)));
    %[Psi,~] = qr(2*rand(n,k)-1,0);
    %isqrt = @(D) diag(sqrt(diag(D)).^-1);
    %Disqrt = isqrt(D);
    %U = Disqrt*Psi;

    %U = rand(n,k);
    %U = ones(n,k);
    % Neumann et al. L1 constraint
    %M = speye(n);
    % Brandt & Hildebrandt
    M = D;

    % Loop over modes
    for i = 1:k
      state = [];
      % TODO: Should probably reformulate as a conic problem
      % 
      % min_U⁺,U⁻   [U⁺' U⁻'] [L -L;-L L] [U⁺; U⁻] + μ(1' M U⁺ + 1' M U⁻)
      % subject to  [c'M  0;0 -c'M ] [U⁺; U⁻] = [1;1]
      %        and  [Uj'M 0;0 -Uj'M] [U⁺; U⁻] = [0;0]  j = 1:n \ i
      %        and  U⁺,U⁻≥ 0
      %

      One = ones(n,1);
      Z = zeros(n,1);
      Uj = U(:,setdiff(1:i-1,i));
      Aeqj = [Uj'*D        -Uj'*D];
      Beqj = [zeros(size(Uj,2),1)];
      fval = inf;
      fvals = [];
      for iter = 1:max_iter

        c = U(:,i)/sqrt(U(:,i)'*D*U(:,i));
        Aeqc = [c'*D   -c'*D];
        Beqc = 1;
        if constrained && i > 1
          % no linear term
          f = zeros(2*n,1);
          Aeq1 = [One'*M One'*M];
          % s := ‖U₁‖₁
          Beq1 = sum(abs(M*U(:,1)));
          Aeq = [Aeq1;Aeqc;Aeqj];
          Beq = [Beq1;Beqc;Beqj];
        else
          f =  mu*[M*One;M*One];
          Aeq = [Aeqc;Aeqj];
          Beq = [Beqc;Beqj];
        end
        Q = [L -L;-L L];
        lx = [Z;Z];
        ux = inf(2*n,1);
        use_quadprog = false;
        if use_quadprog
          params = default_quadprog_param();
          fval_prev = fval;
          [Ui,fval,exitflag,output,lambda] = quadprog( ...
            Q,f, ...
            [],[], ...
            Aeq, Beq, ...
            lx,ux,[],params);
        else
          fval_prev = fval;
          [Ui,state] = ...
            quadprog_box(Q,f,sparse(Aeq),Beq,lx,ux,state,'MaxIter',quadprog_max_iter,'Aeq_li',true);
          fval = 0.5*Ui'*Q*Ui+Ui'*f;
        end
        rel_delta = abs(fval_prev-fval)/abs(fval);
        fvals = [fvals;rel_delta];
        %subplot(2,1,2);
        %semilogy(fvals,'LineWidth',3);
        %drawnow;
        U(:,i) = Ui(1:n)-Ui(n+(1:n));
        if rel_delta < 1e-7
          break;
        end

      end
    end

    Phi = U;
  end
  
  neg = sum(D*Phi)<0;
  % Might as well flip the signs so that each mode is >0 on average
  Phi(:,neg) = -Phi(:,neg);
  % Sort according to primary energy
  [lambda,I] = sort(sum(Phi.*(L*Phi)),'descend');
  Phi = Phi(:,I);
  lambda = diag(lambda);
  
end

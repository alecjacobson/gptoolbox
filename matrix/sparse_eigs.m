function [Phi,lambda] = sparse_eigs(L,D,k,mu)
  % Sparse (aka "compressed") "eigenmodes" minimizing:
  %
  % min_φ  tr( φ' L φ ) + μ ‖φ‖₁  subject to φ' D φ = I
  %
  % This ADMM optimization is implemented according to "Compressed Manifold
  % Modes for Mesh Processing" [Neumann et al. 2014].
  %
  % Inputs:
  %   L  #L by #L sparse, semi-positive definite "stiffness" matrix 
  %   D  #L by #L diagonal positive definite "mass" matrix
  %   k  number of modes requestedd
  %   mu  sparsity parameter (larger is more sparse) {10}
  % Outputs:
  %   Phi  #L by k modes sorted in decreasing order by energy value
  %   lambda  k by k diagonal matrix of energy values (in some way analogous to
  %     eigen values, though L*Phi ≠ lambda*Phi)
  %

  M = D;

  % From Python code.
  tol_abs = 1e-8;
  tol_rel = 1e-6;
  n = size(L,1);
  assert(isdiag(D));
  Dsqrt = diag(sqrt(diag(D)));
  isqrt = @(D) diag(sqrt(diag(D)).^-1);
  Disqrt = isqrt(D);
  % Python codes uses qr factorization of random numbers in [-1,1] range,
  % presummably to start with orthogonal vectors
  %[Psi,~] = qr(2*rand(n,k)-1,0);
  % This is (more) deterministic. The cost is negligible in comparison to the
  % following.
  [Psi,~] = eigs(L,D,k,'sm');

  Phi = Disqrt*Psi;
  S = Phi;
  E = S;
  U = zeros(2*n,k);
  I = speye(n);
  
  rho_prev = nan;
  rho = 1;
  % paper says multiply by n, but python code divides by n.
  % This does seem to be the right units. 
  mu = mu/n;
  
  shrink = @(u,delta) sign(u).*max(abs(u)-delta,0);
  % From Python Code
  check_interval = 10;
  
  num_refactors = 0;
  max_iter = 1000;
  for iter = 1:max_iter
    %if mod(iter,ceil(max_iter/30)) == 2
    %  clf;
    %  pink = [254 194 194]/255;
    %  teal = [144 216 196]/255;
    %  white = [1 1 1];
    %  hold on;
    %  %colormap(isolines_map(cbrewer('RdBu',16)));
    %  colormap(interp1([-1;0;1],[teal;white;pink],linspace(-1,1,16)));
    %  off = 0.7*[max(VV(:,1))-min(VV(:,1)) 0 0];
    %  %order = 1:k;
    %  order = [5 3 4 1 2]; 
    %  t = {};
    %  for ei = 1:numel(order)
    %    R = axisangle2matrix([1 0 0],-pi/2)* ...
    %      axisangle2matrix([0 0 1],-pi/3+(ei-1)*pi/21);
    %    e = order(ei);
    %    t{ei} = tsurf(LF,LV*R+(ei+0.02*ei*ei)*off,'CData',SS*Phi(:,e),fphong,'EdgeColor','none',fsoft);
    %  end
    %  %apply_ambient_occlusion([t{:}],'AddLights',false,'SoftLighting',false,'AO',AO);
    %  camproj('persp');
    %  axis equal;
    %  view(13,15);
    %  %set(gcf,'Color',0.9*[1 1 1]);
    %  set(gcf,'Color',0.4*[1 1 1]);
    %  l = light('Position',[-0.5 10 20],'Style','local');
    %  s = ...
    %    add_shadow([],l,'Fade','local','Color',get(gcf,'Color')*0.9,'BackgroundColor',get(gcf,'Color'));
    %  set(gca,'Visible','off');
    %  set(gca,'Position',[0 0 1 1]);
    %  camlight;
    %  hold off;
    %  add_isolines(t);
    %  figgif('cmm-hand.gif');
    %end
  
    % Update Phi (10)
    UE = U(1:n,:);
    US = U(n+(1:n),:);
    %% Based on the paper:
    %Y = 0.5*(S - US + E - UE);
    %A = Dsqrt*Y;
    %[V,W] = eig(A'*A);
    %% Now  V*W*V' = (Dsqrt*Y)'*(Dsqrt*Y)
    %Phi = Disqrt * (Y*V*isqrt(W)*V');

    %% Based on the Python code:
    %Y = Dsqrt*(E-UE + S-US);
    %% Matlab:  [U,S,V] = svd(A) --> U*S*V' = A
    %% Numpy:  [U,S,VT] = svd(A) --> U*S*VT = A
    %[V,W,VT] = svd(Y'*Y);
    %Phi_prev = Phi;
    %Wisqrt = isqrt(W);
    %Psi = (Y * (V * Wisqrt * VT'));
    %Phi = Disqrt * Psi;
    Psi = rectangular_procrustes(Dsqrt,E-UE + S-US);
    Phi = Disqrt * Psi;
  
    % Python code updates S before E... shouldn't matter
    % Update S (12)
    Sprev = S;
    %% Based on paper (should be same)
    % v = Phi + US;
    % S = sign(v).*max(abs(v) - mu/rho,0);
    % Based on python code
    S = shrink(Phi + US,mu/rho);
  
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
    U = U+[Phi;M*Phi]-[E;S];
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
    % From Python code
    eps_pri = sqrt(k*2)*tol_abs + tol_rel*max(sqrt(sum([Phi;E;S].^2,2)));
    eps_dual = sqrt(k)*tol_abs + tol_rel*rho*sqrt(sum(U(:).^2));
    if residual < eps_pri && dual_residual < eps_dual
      iter
      break;
    end
  end

  neg = sum(D*Phi)<0;
  % Might as well flip the signs so that each mode is >0 on average
  Phi(:,neg) = -Phi(:,neg);
  % Sort according to primary energy
  [lambda,I] = sort(sum(Phi.*(L*Phi)),'descend');
  Phi = Phi(:,I);
  lambda = diag(lambda);
  
end

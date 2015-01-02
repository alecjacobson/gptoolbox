function [U,UT] = isometric_curve_flow(V,varargin)
  % ISOMETRIC_CURVE_FLOW Isometric flow for curves as described by "Robust
  % Fairing via Conformal Curvature Flow" by [Crane et al. 2013]
  %
  % U = isometric_curve_flow(V)
  % U = isometric_curve_flow(V,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   Optional:
  %     'Tao' followed by timestep tao
  % Outputs:
  %   U  #V by 2 list of vertex positions
  %

  % time step
  tao = 1e-2;

  ii = 1;
  while ii < numel(varargin)
    switch varargin{ii}
    case 'Tao'
      assert(ii+1<=numel(varargin));
      ii = ii+1;
      tao = varargin{ii};
    otherwise
      error('Unsupported parameter: %s',varargin{ii});
    end
    ii = ii+1;
  end

  % number of points
  n = size(V,1);

  % Evaluate curvature and edge vectors
  [kappa,alpha,ev,l] = curvature(V);

  % desired flow direction
  kappa_dot = - 2 * kappa;
  
  % defining edge length at each point
  egde_vector = zeros(n,1);
  for k=2:n-1
    edge_vector(k) = 0.5*(norm(V(k+1,:)-V(k,:))+norm(V(k-1,:)-V(k,:)));
  end
  edge_vector(1) = 0.5*(norm(V(2,:)-V(1,:))+norm(V(size(V,1),:)-V(1,:)));
  edge_vector(n) = 0.5*(norm(V(1,:)-V(n,:))+norm(V(n,:)-V(n-1,:)));
  % mass matrix
  M = diag(sparse(0.5*(l([end 1:end-1]) + l))); 
  
  % re-defining orthonormal basis
  C(:,1) = ones(n,1)/sqrt(ones(n,1)'*M*ones(n,1));
  
  C(:,2) = V(:,1) - ((V(:,1)'*M*C(:,1)))*C(:,1);
  C(:,2) = C(:,2)/sqrt(C(:,2)'*M*C(:,2));
  
  C(:,3) = V(:,2) - ((V(:,2)'*M*C(:,1)))*C(:,1) - ((V(:,2)'*M*C(:,2)))*C(:,2) ; 
  C(:,3) = C(:,3)/sqrt(C(:,3)'*M*C(:,3));
  
  % kappa_dot = kappa_dot - ...
  %  (kappa_dot'*M*C(:,1))*C(:,1) - ...
  %  (kappa_dot'*M*C(:,2))*C(:,2) - ...
  %  (kappa_dot'*M*C(:,3))*C(:,3) ;
  kappa_dot = kappa_dot - sum(bsxfun(@times,sum(bsxfun(@times,kappa_dot,M*C)),C),2);

  % Take explicit euler step
  kappa = kappa + tao * kappa_dot;

  % arbitrarily let theta_0 = 0
  theta_0 = 0;
  % % angle between last segment and -x-axis
  %theta_0 = atan2( ...
  %  ev(end,1).*0 - ev(end,2).*-1, ...
  %  ev(end,1).*-1 + ev(end,2).*0 );
  theta_0 = atan2( ...
    -1.*ev(end,2) - 0.*ev(end,1), ...
    -1.*ev(end,1) + 0.*ev(end,2));
  theta = theta_0 + cumsum(kappa .* (0.5 * (l([end 1:end-1]) + l)));

  % recover tangents: T(i,:) is along V(i,:) to V(i+1,:)
  T = bsxfun(@times,l,[cos(theta) sin(theta)]);
  UT = [V(1,:); bsxfun(@plus,V(1,:),cumsum(-T))];

  % Recover positions
  % build laplacian
  L = sparse(1:n,[2:n 1],-1./l,n,n);
  % symmetric
  L = L+L';
  % diagonal diagonals
  L = L - diag(sum(L,2));
  % Bulid rhs
  b = bsxfun(@rdivide,T,l);
  b = b-b([end 1:end-1],:);
  % % solve
  %U = L\b;
  % Fix U(1,:) and solve
  U = V;
  U(2:end,:) = L(2:end,2:end) \ (-L(2:end,1) * U(1,:) + b(2:end,:));

  % % fit a rigid transformation
  %[R,T] = fit_rigid(V,U);
  %U = bsxfun(@plus,U*R,T);

end

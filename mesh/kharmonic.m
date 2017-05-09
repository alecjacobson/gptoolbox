function W = kharmonic(V,F,b,bc,k,varargin)
  % KHARMONIC k-harmonic coordinates, "Harmonic Coordinates for Character
  % Articulation" by Joshi et al, and "An Intuitive Framework for Real-Time
  % Freeform Modeling" section 3 "Precomputed basis functions" by Botsch and
  % Kobbelt, implemented for k>2 using "Mixed Finite Elements for Variational
  % Surface Modeling" [Jacobson et al. 2010]
  %
  % W = kharmonic(V,F,b,bc)
  % W = kharmonic(V,F,b,bc,k);
  % W = kharmonic(V,F,b,bc,k,'ParameterName',ParameterValue);
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %  k  power of laplacian {k=1 %harmonic}
  %    Optional:
  %      'Condensed' followed by whether to use the statically condensed
  %      (positive definite, therefore faster) system rather than the full
  %      KKT system (symmetric though not positive definite, but more stable
  %      numerically).
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions, biharmonic_bounded
  %

  if ~exist('k','var')
    % default to harmonic coordinates
    k = 1;
  end

  condensed = true;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Condensed'}, {'condensed'});
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

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);

  % Build discrete laplacian and mass matrices used by all handles' solves
  if n == 0
    A = adjacency_matrix(F);
    L = A-diag(sum(A,2));
    n = max(F(:));
    M = speye(n);
  else
    L = cotmatrix(V,F);
    M = massmatrix(V,F);
    % NORMALIZE MASSMATRIX (THIS IS IMPORTANT!!)
    M = M./max(abs(diag(M)));
  end

  if condensed
    % build k-laplacian, quadratic coefficients
    switch k
    case 1
      Q = -L;
    case 2
      Q = L*(M\L);
    case 3
      Q = L*((M\L)*(M\L));
      Q = 0.5*(Q+Q');
    otherwise
      Q = -L;
      for ii = 2:k
        Q = Q*(M\-L);
      end
    end

    % Minimize W'QW subject to W(b,:) = bc
    W = min_quad_with_fixed(Q,zeros(n,1),b,bc);
  else
    b = reshape(b,[],1);
    u = reshape(setdiff(1:n,b),[],1);
    a = [b;u];
    Z = 0*M;
    % uses [Jacobson et al. 2010] but exchanges first and second block-columns
    % so that primary variable comes first (also affects solution vector), and
    % then exchange first and second row-blocks to be symmtric (also affects
    % rhs)
    switch k
    % Left commented-out for posterity's sake
    %case 1
    %  % Could use chol but ldl will do...
    %  Q = L(u,u);
    %  B = -L(u,b)*bc;
    %case 2
    %  Q = [ ...
    %     Z(u,u)  L(u,a); ...
    %     L(a,u) -M(a,a); ];
    %  B = ( [Z(u,1:size(bc,2)); -L(a,b) * bc] );
    %case 3
    %  Q = [ ...
    %     Z(u,u)  Z(u,a) L(u,a); ...
    %     Z(a,u)  L(a,a),-M(a,a); ...
    %     L(a,u),-M(a,a) Z(a,a)  ];
    %  B = [Z(u,1:size(bc,2));Z(a,1:size(bc,2));-L(a,b)*bc];
    otherwise
      % Create anti-diagaonal of L's and an anti-off-diagonal of M's
      Q = fliplr(repdiag(fliplr(L(a,a)),k)) + ...
        blkdiag(Z,fliplr(repdiag(fliplr(-M(a,a)),k-1)));
      % Kill first b rows/columns
      Q = Q(numel(b)+1:end,numel(b)+1:end);
      % Zeros followed by -L*bc
      B = [zeros(n*(k-1),size(bc,2));-L(a,b)*bc];
      % Kill first b rows/columns
      B = B(numel(b)+1:end,:);
    end
    % ldl is _much_ faster than backslash
    [s.L,s.D,s.P,s.S] = ldl(Q);
    Wsol = s.S * (s.P * (s.L'\(s.D\(s.L\(s.P' * (s.S * B))))));
    W = zeros(size(V,1),size(bc,2));
    W(b,:) = bc;
    W(u,:) = Wsol(1:numel(u),:);
  end


  %Aeq = sparse(1:numel(b),b,1,numel(b),size(V,1));
  %for c = 1:size(bc,2)
  %  W(:,c) = quadprog(Q,zeros(n,1),[],[],Aeq,bc(:,c));
  %end
end

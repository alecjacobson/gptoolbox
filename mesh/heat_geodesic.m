function [D,u,X,div_X,phi,pre,B,t] = heat_geodesic(varargin)
  % HEAT_GEODESIC  geodesic distance approximation following the method of
  % "Geodesics in Heat" [Krane et al. 2013] D from all points V in the domain
  % (V,F) to a source point/set of points, gamma. The method is motivated by
  % heat using a time parameter t to guide the heat diffusion.
  % 
  % D = heat_geodesic(V,F,gamma,t)
  % [D,u,X,div_X,phi,pre] = heat_geodesic(V,F,gamma,t,'ParamName',ParamValue)
  %
  % Inputs:
  %   V  #V by 3 set of vertex positions 
  %   F  #F by 3 set of face indices
  %     or 
  %      #T by 4 set of tetrahedra indices
  %   gamma  #gamma list of vertex indices of source points
  %   t  time parameter
  %   Optional:
  %     'BoundaryConditions'  followed by one of the following strings
  %       {'average'}: 
  %       'dirichlet'  domain boundaries set to 0 when computing heat
  %         diffusion
  %       'neumann'  heat diffusion solved with implicit neumann conditions
  %       'average'  Uniform average of neumann and Dirichlet solutions
  %       'robin'  The original ArXiv paper called the method of averaging the
  %         Dirichlet and Neumann solutions 'robin'. This was later corrected
  %         (this solution cannot be created by boundary conditions to the heat
  %         equation alone: average solutions is a non-linear operation).  In
  %         any case, 'robin' is treated as 'average'.
  %     'Precomputation' Followed by pre struct as returned by this function
  %     'IntrinsicDelaunay'  followed by whether to use intrinsic Delaunay
  %       Laplacian {false}
  %     'Legacy' followed by bool telling whether to use Alec's legacy
  %       implmentation. In particular this is useful because the original
  %       paper is unclear how the boundary of the domain should be handle with
  %       respect to the seed locations (gamma). What if they overlap? What
  %       should be the boundary conditions for the final poisson solve?
  % Outputs:
  %   D  #V list of geodesic distances from vertices to gamma
  %   u  #V list of results of heat diffusion step
  %   X  #F by 3 list of reversed normalized gradients of u
  %   div_X  #V list od divergence of X
  %   phi  #V list of solution to final poisson equation solve
  %
  D=[];
  % Note: This is Alec's previous implmentation. He is still convinced that
  % this is more exact/correct despite not what's written in the article (and
  % implying that refactoring is necessary).
  legacy = false;
  use_intrinsic = false;

  % mandatory input
  V = varargin{1};
  F = varargin{2};

  % number of domain vertices
  n = size(V,1);
  % simplex size
  ss = size(F,2);
  gamma = varargin{3};
  if nargin>=4 && ~isempty(varargin{4})
    t = varargin{4};
  else
    switch ss
    case 3
      % Section 3.1.1
      AM = sum(doublearea(V,F))/2;
      sF = size(F,1);
      c = 5;
      t = c * AM / sF;
      %t = 20*mean(doublearea(V,F));
    case 4
      AM = sum(volume(V,F));
      sF = size(F,1);
      c = 1e3;
      t = c * AM / sF;
      %t = 20*mean(volume(V,F));
    end
  end

  % option parameter default values
  bc_type = 'average';
  pre = [];
  pre.L = [];
  % precomputation for Dirichlet solve
  pre.D = [];
  % precomputation for Neumann solve
  pre.N = [];
  % precomputation for Poisson solve
  pre.poisson = [];
  u = [];

  ii = 5;
  while(ii <= nargin)
    switch varargin{ii}
    case 'BoundaryConditions'
      ii = ii + 1;
      assert(ii<=nargin);
      bc_type = varargin{ii};
    case 'Precomputation'
      ii = ii + 1;
      assert(ii<=nargin);
      % skip empty input
      if ~isempty(varargin{ii})
        pre = varargin{ii};
      end
    case 'Legacy'
      ii = ii + 1;
      assert(ii<=nargin);
      legacy = varargin{ii};
    case 'u'
      ii = ii + 1;
      assert(ii<=nargin);
      u = varargin{ii};
    case 'IntrinsicDelaunay'
      ii = ii + 1;
      assert(ii<=nargin);
      use_intrinsic = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  % Algorithm 1 in Section 3:
  %   Integrate the heat flow ?? = ??u for time t
  %   Evaluate the vector field X = -???u/|???u|
  %   Solve the Poisson equation ???? = ??????X


  % Integrate the heat flow: "We discretize the heat equation from step I of
  % Algorithm 1 using a single backward Euler step"


  % "... with initial conditions u0 = ??(x)"
  % "Note that a Dirac delta appears as a literal one in this system since we
  % are effectively working with integrated quantities"
  u0 = zeros(n,1);
  u0(gamma) = 1;
  

  % "where ???????? is one third the area of all triangles incident on vertex ...
  % where ???? ??? R|????|??|????| is a diagonal matrix containing the vertex areas"
  if isempty(pre.L)
    if use_intrinsic
    pre.L = intrinsic_delaunay_cotmatrix(V,F);
    else
    pre.L = cotmatrix(V,F);
    end
    pre.M = massmatrix(V,F,'barycentric');
    pre.Q = pre.M - t*pre.L;
    pre.G = grad(V,F);
    pre.Div = div(V,F);
    pre.Tik = 1e-8*speye(size(V,1));
    pre.P = -pre.L*0.5 + pre.Tik;
    % get outline ("boundary") of mesh 
    switch ss
    case 3
      pre.out = unique(reshape(outline(F),[],1));
    case 4
      pre.out = unique(boundary_faces(F));
    end
  end

  % 11/7/2017: Keenan writes that the mass matrix should not be here.
  % B = M*u0;
  B = u0;

  if isempty(u)
    if strcmp(bc_type,'dirichlet') || strcmp(bc_type,'average') || strcmp(bc_type,'robin')

      if legacy
        % remove any gamma from outline
        pre.out = setdiff(pre.out(:),gamma(:));
        % find all boundary vertices and append to gamma
        b = [gamma(:); pre.out]; 
        % "zero Dirichlet conditions"
        bc = [ones(numel(gamma),1); zeros(numel(b)-numel(gamma),1)];
        [uD,pre.D] = min_quad_with_fixed(pre.Q,B,b,bc,[],[],pre.D);
      else
        b = pre.out;
        % Q: How should we deal with gamma ??? out ?
        % This gives *reasonable* results, but doesn't look right for fixing
        % entire boundary. Probably more because of the lack of correct boundary
        % conditions for the final poisson solve.
        bc = -1*ismember(pre.out,gamma);
        %assert(~any(ismember(gamma,out)));
        [uD,pre.D] = min_quad_with_fixed(pre.Q,B,b,bc,[],[],pre.D);
      end
    end
    if strcmp(bc_type,'neumann') || strcmp(bc_type,'average') || strcmp(bc_type,'robin')
      % See Figure 9, pretty sure these are implicit ???x/???n = 0 neumann
      % conditions, though it is not stated in the text. At best, "Neumann
      % conditions prevent heat from flowing out of the domain..."
      if legacy
        b = gamma(:);
        bc = ones(numel(gamma),1);
        [uN,pre.N] = min_quad_with_fixed(pre.Q,B,b,bc,[],[],pre.N);
      else
        [uN,pre.N] = min_quad_with_fixed(pre.Q,B,[],[],[],[],pre.N);
      end
    end

    if strcmp(bc_type,'natural')
      K = keenan(V,F);
      Q = M - t*K;
      B = M*u0;
      G = grad(V,F);
      switch size(F,2)
      case 4
        vol = volume(V,F);
      case 3
        vol = doublearea(V,F);
      end
      vol = vol/sum(vol);
      % kill off affine functions
      A = [
        sum(M)/sum(M(:)); ...
        kron(speye(size(V,2)),vol)'*G];
      uT = min_quad_with_fixed(Q,B,[],[],A,[0;zeros(size(V,2),1)]);
    end

    switch bc_type
    case 'dirichlet'
      u = uD;
    case 'neumann'
      u = uN;
    case {'average','robin'}
      % "We advocate the use of the Robin boundary conditions obtained by taking
      % the mean of the Neumann solution uN and the Dirichlet solution uD, i.e.,
      % u = 0.5*(uN + uD)"
      u = 0.5*(uN+uD);
    case 'natural'
      u = uT;
    otherwise
      error(['Unsupported BoundaryCondtions value: ' bc_type]);
    end
  end

  robust_X = false;
  if robust_X
    X = robust_heat_gradients(V,F,uN,uD,t);
    warning('robust');
  else
    % Evaluate the vector field X
    grad_u = reshape(pre.G*u,size(F,1),size(V,2));
    %grad_u_norm = sqrt(sum(grad_u.^2,2));
    %% normalize grad_u
    %normalized_grad_u = bsxfun(@rdivide,grad_u,grad_u_norm);
    %% correct any zero-norm gradients
    %normalized_grad_u(grad_u_norm == 0,:) = 0;
    % normalize reverse direction; normalizerow will use robust norm that
    % avoids underflow.
    X = -normalizerow(grad_u);
    X(isnan(X)) = 0;
  end
  
  % Solve the Poisson equation 
  % divergence of X
  div_X = pre.Div*X(:);
  if legacy
    [phi,pre.poisson] = min_quad_with_fixed( ...
      -pre.L*0.5,2*div_X,gamma,zeros(numel(gamma),1),[],[],pre.poisson);
  else
    [phi,pre.poisson] = ...
      min_quad_with_fixed(pre.P,2*div_X,[],[],[],[],pre.poisson);
    %Aeq = sparse(ones(1,size(V,1)));
    %Beq = 0;
    %[phi,pre.poisson] = ...
    %  min_quad_with_fixed(-pre.L*0.5,2*div_X,[],[],Aeq,Beq,pre.poisson);
  end
  D = phi;
  % "Note that ???? is unique only up to an additive constant and should be
  % shifted such that the smallest distance value is zero."
  %
  % D should be zero at gamma
  D = D - mean(D(gamma));
  % Flip sign so farthest point is positive distance
  [~,mi] = max(abs(D));
  D = D*sign(D(mi));

end

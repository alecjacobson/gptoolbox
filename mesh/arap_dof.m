function [L,new_C] = arap_dof(varargin)
  % ARAP_DOF compute automatic dof (linear transformation, for now) for point
  % handles in a linear blend skinning deformation. Provided are only the
  % displacements (translations) at each control point.
  % 
  % [L] = arap_dof(V,F,W,C,P,BE,new_C)
  % [L] = arap_dof(V,F,W,C,P,BE,new_C,'PropertyName','PropertyValue',...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of face indices
  %   W  #V by #C list of skinning weights
  %   C  #C by dim list of control vertex positions
  %   P  #P list of indices into C for point controls, { 1:size(C,1) }
  %   BE  #BE by 2 list of bones, pairs of indices into C, connecting control 
  %      vertices in C
  %   %b  #C list of control point indices in V
  %   new_C  #C by dim list of control point pose positions
  %   Optional properties:
  %     'L0'
  %       dim by dim by #C list of linear transformations initial guesses,
  %       optional (default is to use identity transformations)
  %%     'CovarianceScatterMatrix'
  %%       followed by CSM a dim*n by dim*n sparse matrix, see Outputs
  %     'Groups'
  %       followed by #V list of group indices (1 to k) for each vertex, such 
  %       that vertex i is assigned to group G(i)
  %     'Tol'
  %       stopping critera parameter. If variables (linear transformation matrix
  %       entries) change by less than 'tol' the optimization terminates,
  %       default is 0.75 (weak tolerance)
  %     'MaxIter'
  %       max number of local-global iterations, default is 10
  %     'Free'
  %       List of "free" control point, that is, LBS control points that are not controled
  %       by the user. In the optimization this means that interpolation is
  %       released as a constraint for these control points. Indexes C, so to
  %       have a totally free bone you need to specify both joints as free
  %     'Fixed'
  %       List of *totally* fixed handles and their corresponding linear
  %       transformations as a dim by dim by #fixed matrix. Note that this is
  %       not necessarily the complement of 'Free'
  %     'RecenterTranslations' centers translations to control handles before
  %       starting optimization, default is false
  %     'Dynamic' 
  %        #V by dim list of external forces
  %     'TimeStep'
  %        scalar time step value
  %     'Lm1'
  %       dim by dim by #C list of linear transformations at time t-1,
  %       optional (default is to use L0)
  % Outputs:
  %   L dim * dim * #C list of linear transformations, so that new_V = MLBS*L,
  %     where MLBS = lbs_matrix(V,W);
  %%   CSM dim*n by dim*n sparse matrix containing special laplacians along the
  %%     diagonal so that when multiplied by repmat(U,dim,1) gives covariance 
  %%     matrix elements, can be used to speed up next time this function is
  %%     called, see function definitions
  %
  
  % turn off blank constraints warning
  old_warn_blank = warning('QUERY','min_quad_with_fixed:blank_eq');
  warning('OFF','min_quad_with_fixed:blank_eq');
  old_warn_sparse = warning('QUERY','min_quad_with_fixed:sparse_system_dense_constraints');
  warning('OFF','min_quad_with_fixed:sparse_system_dense_constraints');

  % parse input
  V = varargin{1};
  F = varargin{2};
  W = varargin{3};
  C = varargin{4};
  P = varargin{5};
  BE = varargin{6};
  new_C = varargin{7};

  % number of mesh (domain) vertices
  n = size(V,1);
  % dimension of mesh
  dim = size(V,2);
  %warning('Snapping C to V');
  assert(size(C,2) == dim);
  energy = 'spokes';

  % compute distance from every vertex in the mesh to every control vertex
  snap_D = permute(sum((repmat(V,[1,1,size(C,1)]) - ...
    permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
  % use distances to determine closest mesh vertex to each control position
  % such that V(b(i),:) is the closest position to C(i,:) and also W(b(i),:)
  % gives (approximately) the weights near position C(i,:)
  [XXX,b] = min(snap_D);

  % dimension of control points and displacements should match dimension of mesh
  assert(size(C,2) == dim);
  assert(size(new_C,2) == dim);

  % number of point handles
  np = numel(P);
  nb = size(BE,1);

  % weights matrix size should be n by m
  assert(size(W,1) == n && size(W,2) == (np+nb));

  G = [];
  max_iterations = 10;
  tol = 0.75;% tol = 0.75 --> max_iterations ~< 10
  recenter_translations = false;

  free = [];
  fixed = [];
  fixed_L = zeros(dim,dim,0);
  CSM = [];
  % default is no dynamics
  dynamic = false;
  % default is no external forces
  fext = zeros(size(V));
  % defaults is unit time step
  h = 1;
  k = [];


  ii = 8;
  while(ii <= nargin)
    switch varargin{ii}
    case 'L0'
      ii = ii + 1;
      assert(ii<=nargin);
      L0 = varargin{ii};
    %case 'CovarianceScatterMatrix'
    %  ii = ii + 1;
    %  assert(ii<=nargin);
    %  CSM = varargin{ii};
    case 'Groups'
      ii = ii + 1;
      assert(ii<=nargin);
      G = varargin{ii};
      if isempty(G)
        k = [];
      else
        k = max(G);
      end
    case 'Tol'
      ii = ii + 1;
      assert(ii<=nargin);
      tol = varargin{ii};
    case 'MaxIter'
      ii = ii + 1;
      assert(ii<=nargin);
      max_iterations = varargin{ii};
    case 'Free'
      ii = ii + 1;
      assert(ii<=nargin);
      free = varargin{ii};
    case 'Fixed'
      ii = ii + 1;
      assert(ii<=nargin);
      fixed = varargin{ii};
      ii = ii + 1;
      assert(ii<=nargin);
      fixed_L = varargin{ii};
    case 'RecenterTranslations'
      ii = ii + 1;
      assert(ii<=nargin);
      recenter_translations = varargin{ii};
    case 'Dynamic'
      ii = ii + 1;
      dynamic = true;
      assert(ii<=nargin);
      fext = varargin{ii};
    case 'TimeStep'
      ii = ii + 1;
      assert(ii<=nargin);
      h = varargin{ii};
    case 'Lm1'
      ii = ii + 1;
      assert(ii<=nargin);
      Lm1 = varargin{ii};
    case 'Energy'
      ii = ii + 1;
      assert(ii<=nargin);
      energy = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  assert(dynamic || numel(free) < (np+nb));
  assert(numel(fixed) <= (np+nb));
  
  % indices of handles: first points then bones
  indices = 1:(np+nb);

  %% list of indices to handles whose controls are interpolated
  %interpolated = indices(~ismember(indices,free));
  %% need to repeat interpolated indices for each dimension
  %interpolated_dim = [];
  %for d = 0:(dim-1)
  %  interpolated_dim = [interpolated_dim d*m+(interpolated)];
  %end

  %% repeat indices in b for each dimension
  %b_dim = [];
  %for d = 0:(dim-1)
  %  b_dim = [b_dim d*n+(b)];
  %end

  % need to repeat fixed indices for each entry in dim by dim+1 transformation
  % matrix
  fixed_dim = [];
  for d = 0:(((dim+1)*dim)-1)
    fixed_dim = [fixed_dim d*(np+nb)+(fixed)];
  end

  % max group id should be less than or equal to number of groups
  assert(isempty(G) || max(G) <= k);


  % set intial guess
  if(~exist('L0','var') || isempty(L0))
    % stack of identity transformations
    Istack = repmat([eye(dim,dim) zeros(dim,1)],[1 1 (np+nb)]);
    % collect transformations into column
    L = reshape(permute(Istack,[3 1 2]),(np+nb)*dim*(dim+1),1);
  else
    L = L0;
  end

  % handle case where all handles are fixed
  assert(~all(ismember(indices,fixed)));

  % collect fixed transformations into column
  fixed_L = reshape(permute(fixed_L,[3 1 2]),numel(fixed)*dim*(dim+1),1);

  % prefactorization for min_quad_with_fixed
  if(~exist('preF','var'))
    preF = [];
  end
  if(~exist('sym','var'))
    sym = [];
  end


  % force initial guess to interpolate control points, seems like this at best
  % saves one or two iterations
  if recenter_translations
    % doesn't support bones
    assert(nb == 0);
    L = permute(reshape(L,[np dim dim+1]),[2 3 1]);
    L(:,dim+1,:) = (stacktimes(L(:,1:dim,:),-C(P,:)) + new_C(P,:))';
    L = reshape(permute(L,[3 1 2]),np*dim*(dim+1),1);
  end


  % Goal:
  %  Find linear transformations at each control point such that when they're
  %  used in LBS the deformation produced is "as-rigid-as-possible".  In other
  %  words, we want a deformation be where the jacobian of the deformation at
  %  each domain vertex is as close as possible to a rotation.
  %
  %  We want to compute:
  %    argmin ||Jj - Rj||^2 
  %      L
  %  where Jj is the jacobian of the LBS deformation that takes vj --> vj' as
  %  presribed by V' = LBS(W,V,L,D), where L is the set of linear
  %  transformations (the unknowns), D is the set of translations at each
  %  control handle (known). And finally Rj is the closest rotation to Jj.
  %
  %  The first question is how to define this norm, ||Jj - Rj||^2. We can do
  %  this in the same way as Sorkine and Alexa 2007. That is pose the problem
  %  as a linear least squares problem, where we want to minimize the following
  %  energy:
  %
  %  E = ???   ??? ||(vi' - vj') - Ri (vi - vj)||^2    (1)
  %     i???V j???N(i)
  %
  % But notice that vi' and vj' are NOT free variables as they are in Sorkine
  % and Alexa. Instead they are a function of L, namely V' = LBS(W,V,L,D).
  % Because V' is a linear combination of the unknowns L, the complexity of our
  % energy is only reduced.
  %
  % Recall from Sorkine and Alexa, the best rotatitions Ri are well defined
  % functions of V' and thus in our case become a functions of L. 
  % 
  % We may minimize (1) by the same iterative method of Sorkine and Alexa. That
  % is we first treat L as fixed and solve (1) for Ri. Then fixing these Ri
  % solve (1) for L. The continous until convergence is reached. Note as in
  % Sorkine and Alexa each step cannot increase the energy so a converging to a
  % local minimum is guaranteed.
  %
  % So to start the iterations we fix L = L0, perhaps the identity or their
  % values from the last frame. Now V' are known values so we may solve for Ri
  % in (1) in the same manner as Sorkine and Alexa (SVD) or better by polar
  % decomposition.
  %
  % Now that we have values for Ri, call them R0. We solve (1) for new L
  % values. Notice when we fix Ri = R0, (1) becomes quadratic and finding its
  % minimum is straight-forward.
  %
  % First let's write V' = LBS(W,V,L,D) in matrix form, so it's clear V' is
  % just a linear combination of L
  %
  % ...
  % 
  % So know we have V' = MLBS * L, where L is a column vector containing
  % dim*dim*m unknowns.
  %
  % Because both Ri and Vi are known we can write Ri (vi - vj) = bij
  %
  %  E = ???   ??? ||(vi' - vj') - bij||^2    (2)
  %     i???V j???N(i)
  %
  % We can write (1) in matrix form:
  %
  % E = -1 V' * Lapl * V  + 2 V' * Ab * 1 + 1 * Abb * 1  (3)
  % 
  % where Lapl is the discrete mesh laplacian and Ab is an adjacency matrix
  %
  % Abij = cwij * (bji - bij)
  %
  % Abij = cwij * bij * bij
  %
  % Then we may replace V' with MLBS * L:
  %
  % E = -1 (MLBS * L)' * Lapl * (MLBS * L) +
  %     2 (MLBS * L)' * Ab * 1 + 
  %     1 * Abb * 1
  %
  % E = -1 L' * MLBS' * Lapl * MLBS * L +
  %     2 L' * MLBS' * Ab * 1 + 
  %     1 * Abb * 1
  %
  % E = L' * Quad * L + L' * Lin  + C
  %
  % where
  %
  % Quad = -1 * MLBS' * Lapl * MLBS
  % Lin = 2 * MLBS' * Ab * 1
  % C = 1 * Abb * 1
  % 

  % compute matrix MLBS, such that U = reshape(MLBS*L,[n dim]) gives 
  % U = lbs(V,W,L)
  MLBS = lbs_matrix(V,W);

  % cotangent laplacian matrix
  Lcot = cotmatrix(V,F);
  % laplacian matrix
  Lapl = 2 * Lcot;

  %% indices of translations in L, so that L(known) == D
  %known = (dim*dim*m+1):(dim*(dim+1)*m);

  % default 'no-groups' G_sum
  switch energy
  case 'elements'
    G_sum = speye(size(F,1));
    k = size(F,1);
  case {'spokes','spokes-and-rims'}
    G_sum = speye(size(V,1));
    k = size(V,1);
  end
  % fit rotations to each deformed vertex
  if(~exist('CSM','var') || isempty(CSM))
    CSM = covariance_scatter_matrix(V,F,'Energy',energy);
    % if there are groups then condense scatter matrix to only build
    % covariance matrices for each group
    if ~isempty(G)
      G_sum = group_sum_matrix(G,k);
      %CSM = [G_sum sparse(k,n); sparse(k,n) G_sum] * CSM;
      CSM = repdiag(G_sum,dim)  * CSM;
    end
  end

  % precompute CSM times MLBS for each dimension
  CSM_MLBS = cell(dim,1);
  for ii = 1:dim
    %CSM_MLBS{ii} = CSM*repmat(MLBS((ii-1)*n+(1:n),:),dim,1);
    CSM_MLBS{ii} = sparse(size(CSM,1),size(MLBS,2));
    Mi = MLBS((ii-1)*n+(1:n),:);
    for jj = 1:dim
      CSM_MLBS{ii}((jj-1)*k+(1:k),:) = ...
        CSM((jj-1)*k+(1:k),(jj-1)*n+(1:n)) * ...
        Mi;
    end
  end


  % precompute arap_rhs matrix
  [XXX,K] = arap_rhs(V,F,[],'Energy',energy);
  MLBS_KG = -4 * MLBS' * K * repdiag(G_sum',dim*dim);
  
  % precompute system matrix
  A = repdiag(-Lapl,dim);
  AA = MLBS'*A*MLBS;
  if(any(any(AA-AA')))
    oldAA = AA;
    % "fix" numerical non symmetry
    AA = 0.5*(AA+AA');
    %max(abs(AA-AA'))
    %max(abs(AA-oldAA))
    %error
  end

  if dynamic
    L0 = L;
    if ~exist('Lm1','var')
      Lm1 = L0;
    end
    V0 = MLBS * L0;
    Vm1 = MLBS * Lm1;
    M = massmatrix(V,F,'voronoi');
    DQ = MLBS'*(0.5*1/h^3*repdiag(M,dim))*MLBS;
    Dl = MLBS'*(1/h^3*repdiag(M,dim)*(-2*V0(:) + Vm1(:)) + fext(:));
  end

  iteration = 0;
  % loop until reached max number of iterations or change in solution is less
  % than stopping parameter
  while( iteration < max_iterations && ...
    (iteration == 0 || max(abs(L(:)-L_prev(:)))>tol))
    % keep track of previous solution
    L_prev = L;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Local step: fix positions, fit rotations per mesh vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% rearrange L from a column to a stack of transformations
    S = zeros(k*dim,dim);
    for ii = 1:dim
      S(:,ii) = CSM_MLBS{ii}*L;
    end
    % dim by dim by n list of covariance matrices
    S = permute(reshape(S,[k dim dim]),[2 3 1]);
    %U = reshape(MLBS*L,[n dim]);
    % compute covariance matrix elements
    R = fit_rotations(S);

    % This is still the SLOW way of building the right hand side, the entire
    % step of building the right hand side should collapse into a
    % #handles*dim*dim+1 by #groups matrix times the rotations at for group
    
    %% distribute group rotations to vertices in each group
    %if ~isempty(G)
    %  R = R(:,:,G);
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % "Global" step: fix rotations per mesh vertex, solve for
    % linear transformations at handles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %rhs_B = arap_rhs(V,F,R,Lcot);

    % E = ???k:dim -V'*Lapl*V + V'*-4*rhs_B + C
    % E = V(:)'*A*V(:) + V(:)'*B + C
    %precomputed A = repdiag(-Lapl,dim);
    %B = -4*rhs_B(:);
    BB = MLBS_KG * reshape(permute(R,[3 1 2]),k*dim*dim,1);
    % Now replace V(:) with MLBS * L
    % E = (MLBS*L)'*A*(MLBS*L) + (MLBS*L)'*B + C
    % E = L'*MLBS'*A*MLBS*L + L'*MLBS'*B + C
    % E = L'*AA*L + L'*BB + C
    %precomputed AA = MLBS'*A*MLBS;
    %BB = MLBS'*B;
    % ???E/???L = 2*AA*L + BB
    % Solve with L = (2*AA)\-BB
    % L = (2*AA)\-BB;
    % BUT we need to add constraints
    %L = min_quad_with_fixed(AA,BB,known,D(:));

    %% BUT we need to ensure that control points are interpolated
    %% linear equality constraint constant values: each coordinate of each
    %% control point's rest position
    %Beq = C+D;
    %Beq = Beq(:);

    %%% This assumes that to interpolate the dragged position of handle i (at
    %%% point ci) the weight of handle i is 1: wi(ci) = 1 and all other weights
    %%% are 0 there wj(ci) = 0, j != i
    %AeqI = repmat(1:dim*m,1,dim+1)';
    %AeqJ = 1:(m*dim*(dim+1))';
    %AeqV = [reshape(repmat(C,dim,1),m*dim*dim,1);ones(m*dim,1)];
    %Aeq = sparse(AeqI,AeqJ,AeqV,dim*m,dim*(dim+1)*m);
    %% This assumes that all weights are well defined at the positions of
    %% interpolated handles, it is then a generalization of the above version
    %% find mesh points
    %Aeq = MLBS(b_dim,:);

    %% throw away constraints for "free" handles, only keep interpolated
    %% handles' constraints
    %Aeq = Aeq(interpolated_dim,:);
    %Beq = Beq(interpolated_dim,:);

    % total number of linear equality constraints, also the current row in the
    % constraints matrix
    Beq = [];
    Aeq = [];
    % whether weights at a given handle (or along a given bone) can be assumed
    % to be exactly 1 for that handle and exactly for all others, including
    % extras
    weights_perfectly_interpolate_handles = true;
    if weights_perfectly_interpolate_handles
      % constraint for each point handle
      for ii = 1:numel(P)
        if ~ismember(P(ii),free)
          % ii is index in weights
          wi = ii;
          % build interpolating weights for this handle
          w = full(sparse(1,wi,1,1,size(W,2)));
          % vertex position of this handle
          v = C(P(ii),:);
          % build phony LBS matrix for only this handle's weights at this
          % handle
          mlbs = lbs_matrix(v,w);
          % constraint for each dimension
          for d = 1:dim
            Aeq = [Aeq; mlbs(d,:)];
            Beq = [Beq;new_C(P(ii),d)];
          end
        end
      end
      % constraint for each bone
      for ii = 1:size(BE,1)
        % offset weight index from points
        wi = ii + numel(P);
        % build interpolating weights for this handle
        w = full(sparse(1,wi,1,1,size(W,2)));
        % for either endpoint
        for bi = BE(ii,:)
          if ~ismember(bi,free)
            % vertex position of this handle
            v = C(bi,:);
            % build phony LBS matrix for only this handle's weights at this
            % handle
            mlbs = lbs_matrix(v,w);
            % constraint for each dimension
            for d = 1:dim
              Aeq = [Aeq; mlbs(d,:)];
              Beq = [Beq;new_C(bi,d)];
            end
          end
        end
      end
    else
      % add constraint for each point handle
      for pi = P
        if ~ismember(pi,free)
          % for each dimension
          for d = 0:(dim-1)
            % index of this coordinate of pi's position in V/W, pi is an element of
            % P which indexes C
            ci = d*n + b(pi);
            Aeq = [Aeq; MLBS(ci,:)];
            Beq = [Beq;new_C(pi,d+1)];
          end
        end
      end
      % add constraint for each unique bone endpoint
      for bi = unique(BE(:))'
        if ~ismember(bi,free)
          % for each dimension
          for d = 0:(dim-1)
            % start point
            ci = d*n + b(bi);
            Aeq = [Aeq; MLBS(ci,:)];
            Beq = [Beq;new_C(bi,d+1)];
          end
        end
      end
    end

    % regularization
    regularize = false;
    if regularize
      if issparse(AA)
        free_mask = sparse(1,(np+nb));
      else
        free_mask = zeros(1,(np+nb));
      end
      free_mask(free) = 1;
      % don't regularize interpolated or fixed DOF
      Gamma = repdiag(diag(free_mask),(dim+1)*dim);
      AA = AA + 0.1*Gamma;
    end

    % I don't think the following is true nor needed. I think this was just
    % hiding the robustness problem with using a badly conditioned system in
    % min_quad_with_fixed > lu_lagrange, which hopefully now is fixed
    %% AA somehow gets -- seemingly random -- small crap in it during
    %% construction, this is a HACK to get rid of it
    %r = min(abs(AA(find(AA))))*1e-10;
    %AA = round(AA/r)*r;

    if dynamic
      AA = AA + DQ;
      BB = BB + Dl;
    end
    [L,preF] = min_quad_with_fixed(AA,BB,fixed_dim,fixed_L,Aeq,Beq,preF);
    %% uncomment to see transformations
    %permute(reshape(L,[m dim dim+1]),[2 3 1])

    % increment number of iterations
    iteration = iteration + 1;
  end

  % project linear transformations to rotations
  project_to_rotations = false;

  %% convert to stack
  %L = permute(reshape(L,[m dim dim+1]),[2 3 1]);
  %for ii = interpolated
  %  if(project_to_rotations)
  %    [su,ss,sv]=svd(L(:,1:dim,ii)');
  %    L(:,1:dim,ii) = sv*su';
  %  end
  %  % reset translations so control points are interpolated
  %  L(:,dim+1,ii) = L(:,1:dim,ii)*(-C(ii,:)') + C(ii,:)' + D(ii,:)';
  %end
  %% convert to column
  %L = reshape(permute(L,[3 1 2]),m*dim*(dim+1),1);

  %R = rand(dim,dim,n);
  %B = zeros(n,n,dim);
  %%Computes energy in non-vectorized way. For double checking terms
  %L = cotmatrix(V,F);
  %E = 0;
  %for i = 1:n
  %for j = 1:n
  %  if( i ~= j )
  %    B(i,j,:) = R(:,:,i) * ((V(i,:) - V(j,:))');
  %    for k = 1:dim
  %      E = E + L(i,j).*(V(i,k) - V(j,k) -B(i,j,k)).^2;
  %    end
  %  end 
  %end
  %end
  %E


  %AdjCot = (Lcot-diag(diag(Lcot)));
  %E = 0;
  %for k = 1:dim
  %  % diagonal already canceled out
  %  Ab = 2*AdjCot.*(B(:,:,k)'-B(:,:,k));
  %  Abb = AdjCot.*B(:,:,k).*B(:,:,k);
  %  %QE = V'*(-2*Lcot)*V;
  %  %LE = V'*Ab*ones(size(V));
  %  %CE = ones(size(V))'*Abb*ones(size(V));
  %  E = E + ...
  %    -V(:,k)'*Lapl*V(:,k) + ...
  %    V(:,k)'*Ab*ones(n,1) + ...
  %    ones(n,1)'*Abb*ones(n,1);
  %end
  %E

  %E = 0;
  %rhs_B = arap_rhs(V,F,R,Lcot);
  %for k = 1:dim
  %  Abb =   AdjCot.*B(:,:,k).*B(:,:,k);
  %  E = E + ...
  %    -V(:,k)'*Lapl*V(:,k) + ...
  %    V(:,k)'*-4*rhs_B(:,k) + ...
  %    ones(n,1)'*Abb*ones(n,1);
  %end
  %E

  %A = -Lapl;
  %BB = -4*rhs_B;
  %E = 0;
  %LE = 0;
  %for k = 1:dim
  %  Abb = AdjCot.*B(:,:,k).*B(:,:,k);
  %  C = ones(n,1)'*Abb*ones(n,1);
  %  E = E + V(:,k)'*A*V(:,k) + V(:,k)'*BB(:,k) + C;
  %end
  %E


  %AA = repdiag(-Lapl,dim);
  %%%[Ai,Aj,Av] = find(A);
  %%%AAi = Ai;
  %%%AAj = Aj;
  %%%AAv = Av;
  %%%for k = 2:dim
  %%%  AAi = [AAi (k-1)*n+Ai];
  %%%  AAj = [AAj (k-1)*n+Aj];
  %%%  AAv = [AAv Av];
  %%%end
  %%%AA = sparse(AAi,AAj,AAv,n*dim,n*dim);
  %%AA = sparse(n*dim,n*dim);
  %%for k = 1:dim
  %%  d = (k - 1)*n + 1;
  %%  AA(d:(d + n-1), d:(d + n-1)) = A;
  %%end

  %BB = -4*rhs_B(:);
  %C = 0;
  %for k = 1:dim
  %  Abb = AdjCot.*B(:,:,k).*B(:,:,k);
  %  C = C + ones(n,1)'*Abb*ones(n,1);
  %end
  %E = V(:)'*AA*V(:) + V(:)'*BB + C;
  %E


  %%E

  % set free control point's positions in new_C
  for pi = P
    if ismember(pi,free)
      for d = 0:(dim-1)
        new_C(pi,d+1) = MLBS(b(pi)+n*d,:)*L;
      end
    end
  end
  for bi = unique(BE(:))'
    if ismember(bi,free)
      for d = 0:(dim-1)
        new_C(bi,d+1) = MLBS(b(bi)+n*d,:)*L;
      end
    end
  end

  % reset min_quad_with_fixed:blank_eq warning
  warning(old_warn_blank.state,'min_quad_with_fixed:blank_eq');
  warning(old_warn_sparse.state,'min_quad_with_fixed:sparse_system_dense_constraints');
end

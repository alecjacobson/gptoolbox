function [U,U1,UF] = takeo_arap(varargin)
  % TAKEO_ARAP Solve "As-Rigid-As-Possible" according to "As-Rigid-As-Possible
  % Shape Manipulation" by Igarashi et al. Given a rest mesh (V,F) and list of
  % constraint vertex indices (b) and their new postions (bc) solve for pose
  % mesh positions (U)
  %
  % U = takeo_arap(V,F,b,bc) Original implementation of takeo's arap: i.e. use
  %   Least Squares Conformal Map as initial guess followed by "1 iteration" of
  %   rotation fitting and poisson solve.
  %
  % U = takeo_arap(V,F,b,bc,'ParameterName','ParameterValue');
  %
  %
  % Inputs:
  %   V  #V by 2 list of rest domain positions
  %   F  #F by 3 list of triangle indices into V
  %   b  #b list of indices of constraint (boundary) vertices
  %   bc  #b by 2 list of constraint positions for b
  %   Optional:
  %     'V0' #V by dim list of initial guess positions
  %       dim by dim by #C list of linear transformations initial guesses,
  %       optional (default is to use identity transformations)
  %     'Tol'
  %       stopping critera parameter. If variables (linear transformation matrix
  %       entries) change by less than 'tol' the optimization terminates,
  %       default is 0.75 (weak tolerance)
  %     'MaxIter'
  %       max number of local-global iterations, default is 1, so that calling
  %       without specifying iterations is equivalent to takeos original
  %       implementation
  % Outputs:
  %   U  #V by 2 list of new positions
  %   U1  #V by 2 list of positions before final fit and poisson iteration
  %   UF  #F*3 by 2 list of positions after final fit but before final poisson
  %     iteration
  %
  % See also: arap, takeo_asap
  %

  % parse input
  V = varargin{1};
  F = varargin{2};
  b = varargin{3};
  bc = varargin{4};

  % default values
  max_iterations = 1;
  tol = 0.001;

  % average edge length, helps later to keep tolerance somewhat domain
  % independent
  h = avgedge(V,F);

  % number of vertices
  n = size(V,1);
  % number of dimensions
  dim = size(V,2);
  % only works in 2D
  assert(dim == 2)
  % number of triangles
  nt = size(F,1);

  ii = 5;
  while(ii <= nargin)
    if strcmp(varargin{ii},'V0');
      ii = ii + 1;
      assert(ii<=nargin);
      U = varargin{ii};
    elseif strcmp(varargin{ii},'Tol');
      ii = ii + 1;
      assert(ii<=nargin);
      tol = varargin{ii};
    elseif strcmp(varargin{ii},'MaxIter');
      ii = ii + 1;
      assert(ii<=nargin);
      max_iterations = varargin{ii};
    end
    ii = ii + 1;
  end

  if ~exist('U','var') || isempty(U)
    % First step: solve ASAP deformation
    [U,XI,YI] = takeo_asap(V,F,b,bc);
  else
    [XI,YI] = relative_coordinates(V,F);
    % reshape coordinate into single tall column
    XI = reshape(XI,size(F,1)*3,1);
    YI = reshape(YI,size(F,1)*3,1);
  end

  assert(all(size(V) == size(U)));


  % indices of each triangle
  t = 1:nt;
  % indices of three corners of each triangle as seen by other corners
  ti = [t+(1-1)*nt t+(2-1)*nt t+(3-1)*nt];
  tj = [t+(2-1)*nt t+(3-1)*nt t+(1-1)*nt];
  tk = [t+(3-1)*nt t+(1-1)*nt t+(2-1)*nt];
  % indices to each triangle
  i = [F(t,1);F(t,2);F(t,3)]';
  j = [F(t,2);F(t,3);F(t,1)]';
  k = [F(t,3);F(t,1);F(t,2)]';
  % Rename for convenience
  xi = XI(ti);
  yi = YI(ti);
  xj = XI(tj);
  yj = YI(tj);
  xk = XI(tk);
  yk = YI(tk);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Prepare fitting system
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % quadratic part
  QI = [[ti ti+3*nt]'; [ti ti+3*nt]'; [tj tj+3*nt]';[ti ti+3*nt]';[tj+3*nt tj]'];  
  QJ = [[ti ti+3*nt]'; [tj tj+3*nt]'; [ti ti+3*nt]';[tj+3*nt tj]';[ti ti+3*nt]'];  
  QV = [ ...
    (1- 2*[xk;xk] + [xk;xk].^2 + [xj;xj].^2 + [-yk;yk].^2 + [-yj;yj].^2) ;...
    (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...
    (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...  
    ([yk;-yk]); ...
    ([yk;-yk]); ...
    ];
  Qfit = sparse(QI,QJ,QV,2*3*nt,2*3*nt);
  % prepare fitting output
  UF = zeros(3*nt,dim);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Prepare poisson solve
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % cotan weights
  Lcot = cotmatrix(V,F);
  % quadratic coefficients
  QI = [j j i i]';
  QJ = [j i j i]';
  QV = [ ...
       Lcot(sub2ind(size(Lcot),i,j)'); ...
      -Lcot(sub2ind(size(Lcot),i,j)'); ...
      -Lcot(sub2ind(size(Lcot),i,j)'); ...
       Lcot(sub2ind(size(Lcot),i,j)')];
  Qpoi = sparse(QI,QJ,QV,n,n);
  %Q = -Lcot;

  
  iteration = 0;
  U1 = V;
  while( iteration < max_iterations && (iteration == 0 || max(abs(U(:)-U1(:)))>tol*h))
    U1 = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Construct and solve fitting problem
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % repeat each vertex for each triangle corner
    U1F = [U1(F(:,1),:);U1(F(:,2),:);U1(F(:,3),:)];
    % linear part for fitting
    Lfit = ...
      -2*(U1F([tk tk+3*nt]').*(1 - [xk;xk]) + ...
      U1F([tj tj+3*nt]').*([xj;xj]) + ...
      U1F([tk+3*nt tk]').*(-[-yk;yk]) + ...
      U1F([tj+3*nt tj]').*([-yj;yj]));
    % solve
    UF(:) = min_quad_with_fixed(Qfit,Lfit,[],[]);
    % ratio of old size to new size
    ratioF = ...
      sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2))./ ...
      sqrt(sum((UF(t,:)-UF(t+nt,:)).^2,2));
    % repeat for each coordinate
    ratioF = repmat(ratioF,1,2);
    % centers of mass
    cmF = (UF(t,:) + UF(t+nt,:) + UF(t+2*nt,:))/3;
    % scale each triangle to original size around its center of mass
    UF(t,:)      = cmF + ratioF.*(UF(t     ,:) - cmF);
    UF(t+nt,:)   = cmF + ratioF.*(UF(t+nt  ,:) - cmF);
    UF(t+2*nt,:) = cmF + ratioF.*(UF(t+2*nt,:) - cmF);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Construct and solve poisson problem
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % linear coefficients of poison problem
    LI = [i i];
    LJ = [ones(3*nt,1) 2*ones(3*nt,1)];
    LV = [ ...
      2*(repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
        (UF(tj,:) - UF(ti,:)) - ...
      repmat(Lcot(sub2ind(size(Lcot),k,i)'),1,2).* ...
      (UF(ti,:) - UF(tk,:)))];
    Lpoi = sparse(LI(:),LJ(:),LV(:),n,2);
    % solve in each coordinate
    U(:,1) = min_quad_with_fixed(Qpoi,Lpoi(:,1),b,bc(:,1));
    U(:,2) = min_quad_with_fixed(Qpoi,Lpoi(:,2),b,bc(:,2));
    
    % increment iteration counter
    iteration = iteration + 1;
  end

end 

  %U = V0;
  %E = 0;
  %Eq = 0;
  %El = 0;
  %Elxy = 0;
  %Ec = 0;
  %for t = 1:nt
  %  for i = 1:3
  %    j = mod(i,3)+1;
  %    k = mod(i+1,3)+1;
  %    xi = XI(t+nt*(i-1));
  %    yi = YI(t+nt*(i-1));
  %    E = E + ...
  %      sum((U(F(t,j),:) + [xi xi].*(U(F(t,k),:)-U(F(t,j),:)) + [yi -yi].*fliplr(U(F(t,k),:)-U(F(t,j),:)) - V0(F(t,i),:)).^2,2);
  %    Eq = Eq + sum((U(F(t,j),:) + [xi xi].*(U(F(t,k),:)-U(F(t,j),:)) + [yi -yi].*fliplr(U(F(t,k),:)-U(F(t,j),:))).^2,2);
  %    El = El + sum(-2*(U(F(t,j),:) + [xi xi].*(U(F(t,k),:)-U(F(t,j),:)) + [yi -yi].*fliplr(U(F(t,k),:)-U(F(t,j),:))).*(V0(F(t,i),:)),2);
  %    %Elxy = Elxy + sum(-2*(U(F(t,j),1) + [xi].*(U(F(t,k),1)-U(F(t,j),1)) + [yi].*fliplr(U(F(t,k),2)-U(F(t,j),2))).*(V0(F(t,i),1)),2);
  %    %Elxy = Elxy + sum(-2*(U(F(t,j),2) + [xi].*(U(F(t,k),2)-U(F(t,j),2)) + [-yi].*fliplr(U(F(t,k),1)-U(F(t,j),1))).*(V0(F(t,i),2)),2);
  %    %El = El + sum(sum( ...
  %    %  -2*(UF(ti,:) -[xk xk].*(UF(ti,:)) -[yk -yk].*fliplr(UF(ti,:))).*V0F(tk,:) + -2*([xj xj].*UF(ti,:)+ [yj -yj].*fliplr(UF(ti,:))).*(V0F(tj,:)),2));
  %    Ec = Ec + sum(V0(F(t,i),:).^2,2);
  %  end
  %end
  %[Eq El Ec E]

  %UF = [U(F(:,1),:);U(F(:,2),:);U(F(:,3),:)];
  %new_UF = zeros(size(UF));

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %for t = 1:nt
  %  ti = t+([1 2 3]-1)*nt;
  %  tj = t+([2 3 1]-1)*nt;
  %  tk = t+([3 1 2]-1)*nt;
  %  assert(~any(new_UF([ti ti+3*nt])));
  %  xi = XI(ti);
  %  yi = YI(ti);
  %  xj = XI(tj);
  %  yj = YI(tj);
  %  xk = XI(tk);
  %  yk = YI(tk);
  %  % quadratic part
  %  QI = [[ti ti+3*nt]'; [ti ti+3*nt]'; [tj tj+3*nt]';[ti ti+3*nt]';[tj+3*nt tj]'];  
  %  QJ = [[ti ti+3*nt]'; [tj tj+3*nt]'; [ti ti+3*nt]';[tj+3*nt tj]';[ti ti+3*nt]'];  
  %  QV = [ ...
  %    (1- 2*[xk;xk] + [xk;xk].^2 + [xj;xj].^2 + [-yk;yk].^2 + [-yj;yj].^2) ;...
  %    (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...
  %    (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...  
  %    ([yk;-yk]); ...
  %    ([yk;-yk]); ...
  %    ];
  %  Q = sparse(QI,QJ,QV,2*3*nt,2*3*nt);
  %  % linear part
  %  L = -2*(V0F([tk tk+3*nt]').*(ones(3*2,1) - [xk;xk]) + V0F([tj tj+3*nt]').*([xj;xj]) + V0F([tk+3*nt tk]').*(-[-yk;yk])+V0F([tj+3*nt tj]').*([-yj;yj]));
  %  Eq = Eq + UF(:)'*Q*UF(:);
  %  El = El + UF([ti ti+3*nt]')'*L;
  %  % constant part
  %  Ec = Ec + sum(sum(V0F(ti,:).^2,2));
  %  E = E + Eq + El + Ec;
  %end
  %[Eq El Ec E]

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %for t = 1:nt
  %  ti = t+([1 2 3]-1)*nt;
  %  tj = t+([2 3 1]-1)*nt;
  %  tk = t+([3 1 2]-1)*nt;
  %  assert(~any(new_UF([ti ti+3*nt])));
  %  xi = XI(ti);
  %  yi = YI(ti);
  %  xj = XI(tj);
  %  yj = YI(tj);
  %  xk = XI(tk);
  %  yk = YI(tk);
  %  % quadratic part
  %  QI = [[ti ti+3*nt]'; [ti ti+3*nt]'; [tj tj+3*nt]';[ti ti+3*nt]';[tj+3*nt tj]'];  
  %  QJ = [[ti ti+3*nt]'; [tj tj+3*nt]'; [ti ti+3*nt]';[tj+3*nt tj]';[ti ti+3*nt]'];  
  %  QV = [ ...
  %    (1- 2*[xk;xk] + [xk;xk].^2 + [xj;xj].^2 + [-yk;yk].^2 + [-yj;yj].^2) ;...
  %    (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...
  %    (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...  
  %    ([yk;-yk]); ...
  %    ([yk;-yk]); ...
  %    ];
  %  Q = sparse(QI,QJ,QV,2*3*nt,2*3*nt);
  %  % linear part
  %  L = -2*(V0F([tk tk+3*nt]').*(ones(3*2,1) - [xk;xk]) + V0F([tj tj+3*nt]').*([xj;xj]) + V0F([tk+3*nt tk]').*(-[-yk;yk])+V0F([tj+3*nt tj]').*([-yj;yj]));
  %  Eq = Eq + UF([ti ti+3*nt]')'*Q([ti ti+3*nt],[ti ti+3*nt])*UF([ti ti+3*nt]');
  %  El = El + UF([ti ti+3*nt]')'*L;
  %  % constant part
  %  Ec = Ec + sum(sum(V0F(ti,:).^2,2));
  %  E = E + Eq + El + Ec;
  %  new_UF([ti ti+3*nt]') = min_quad_with_fixed(Q([ti ti+3*nt],[ti ti+3*nt]),L,[],[]);
  %end
  %[Eq El Ec E]
  %max(abs([V0F(:)-new_UF(:)]))

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %t = 1:nt;
  %ti = [t+(1-1)*nt t+(2-1)*nt t+(3-1)*nt];
  %tj = [t+(2-1)*nt t+(3-1)*nt t+(1-1)*nt];
  %tk = [t+(3-1)*nt t+(1-1)*nt t+(2-1)*nt];
  %xi = XI(ti);
  %yi = YI(ti);
  %xj = XI(tj);
  %yj = YI(tj);
  %xk = XI(tk);
  %yk = YI(tk);

  %% quadratic part
  %QI = [[ti ti+3*nt]'; [ti ti+3*nt]'; [tj tj+3*nt]';[ti ti+3*nt]';[tj+3*nt tj]'];  
  %QJ = [[ti ti+3*nt]'; [tj tj+3*nt]'; [ti ti+3*nt]';[tj+3*nt tj]';[ti ti+3*nt]'];  
  %QV = [ ...
  %  (1- 2*[xk;xk] + [xk;xk].^2 + [xj;xj].^2 + [-yk;yk].^2 + [-yj;yj].^2) ;...
  %  (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...
  %  (-[xk;xk].^2 +[xk;xk]-[-yk;yk].^2); ...  
  %  ([yk;-yk]); ...
  %  ([yk;-yk]); ...
  %  ];
  %Q = sparse(QI,QJ,QV,2*3*nt,2*3*nt);
  %% linear part
  %L = -2*(V0F([tk tk+3*nt]').*(1 - [xk;xk]) + V0F([tj tj+3*nt]').*([xj;xj]) + V0F([tk+3*nt tk]').*(-[-yk;yk])+V0F([tj+3*nt tj]').*([-yj;yj]));
  %Eq = Eq + UF([ti ti+3*nt]')'*Q([ti ti+3*nt],[ti ti+3*nt])*UF([ti ti+3*nt]');
  %El = El + UF([ti ti+3*nt]')'*L;
  %% constant part
  %Ec = Ec + sum(sum(V0F(ti,:).^2,2));
  %E = E + Eq + El + Ec;
  %new_UF([ti ti+3*nt]') = min_quad_with_fixed(Q([ti ti+3*nt],[ti ti+3*nt]),L,[],[]);
  %[Eq El Ec E]
  %max(abs([V0F(:)-new_UF(:)]))

  %U = V;
  %U = rand(size(V));
  %Lcot = cotmatrix(V,F);
  %E = 0;
  %for t = 1:nt
  %  for i = 1:3
  %    j = mod(i,3)+1;
  %    E = E + ...
  %      Lcot(F(t,i),F(t,j)) * ...
  %      sum(...
  %        (U(F(t,j),:) - U(F(t,i),:) - (UF(t+(j-1)*nt,:) - UF(t+(i-1)*nt,:))).^2 ...
  %     ,2);;
  %  end
  %end
  %E

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %for t = 1:nt
  %  for i = 1:3
  %    j = mod(i,3)+1;
  %    Eq = Eq + ...
  %      sum(...
  %        Lcot(F(t,i),F(t,j)) * ...
  %        (U(F(t,j),:) - U(F(t,i),:)).^2 ...
  %     ,2);
  %    El = El + ...
  %      sum(...
  %        Lcot(F(t,i),F(t,j)) * ...
  %        -2*(U(F(t,j),:) - U(F(t,i),:)).*(UF(t+(j-1)*nt,:) - UF(t+(i-1)*nt,:)) ...
  %     ,2);
  %    Ec = Ec + ...
  %      sum(...
  %        Lcot(F(t,i),F(t,j)) * ...
  %        (UF(t+(j-1)*nt,:) - UF(t+(i-1)*nt,:)).^2 ...
  %     ,2);
  %  end
  %end
  %E = E + Eq + El + Ec;
  %[Eq El Ec E]

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %for t = 1:nt
  %  ti = t+([1 2 3]-1)*nt;
  %  tj = t+([2 3 1]-1)*nt;
  %  i = F(t,[1 2 3]);
  %  j = F(t,[2 3 1]);
  %  Eq = Eq + ...
  %    sum(sum(...
  %      repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %      (U(j,:) - U(i,:)).^2 ...
  %   ,2));
  %  El = El + ...
  %    sum(sum(...
  %    repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %      -2.*(UF(tj,:) - UF(ti,:)).* ...
  %      (U(j,:) - U(i,:)) ...
  %   ,2));
  %  Ec = Ec + ...
  %    sum(sum(...
  %    repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %      (UF(tj,:) - UF(ti,:)).^2 ...
  %    ,2));
  %end
  %E = E + Eq + El + Ec;
  %[Eq El Ec E]

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %t = 1:nt;
  %ti = [t+(1-1)*nt t+(2-1)*nt t+(3-1)*nt];
  %tj = [t+(2-1)*nt t+(3-1)*nt t+(1-1)*nt];
  %i = [F(t,1);F(t,2);F(t,3)]';
  %j = [F(t,2);F(t,3);F(t,1)]';
  %Eq = Eq + ...
  %  sum(sum(...
  %    repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %    (U(j,:) - U(i,:)).^2 ...
  % ,2));
  %El = El + ...
  %  sum(sum(...
  %  repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %    -2.*(UF(tj,:) - UF(ti,:)).* ...
  %    (U(j,:) - U(i,:)) ...
  % ,2));
  %Ec = Ec + ...
  %  sum(sum(...
  %  repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %    (UF(tj,:) - UF(ti,:)).^2 ...
  %  ,2));
  %E = E + Eq + El + Ec;
  %[Eq El Ec E]

  %E = 0;
  %Eq = 0;
  %El = 0;
  %Ec = 0;
  %t = 1:nt;
  %ti = [t+(1-1)*nt t+(2-1)*nt t+(3-1)*nt];
  %tj = [t+(2-1)*nt t+(3-1)*nt t+(1-1)*nt];
  %tk = [t+(3-1)*nt t+(1-1)*nt t+(2-1)*nt];
  %i = [F(t,1);F(t,2);F(t,3)]';
  %j = [F(t,2);F(t,3);F(t,1)]';
  %k = [F(t,3);F(t,1);F(t,2)]';
  %Eq = Eq + ...
  %  sum(sum(...
  %    repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).* ...
  %    (U(j,:).^2 + -2*U(j,:).*U(i,:)+ U(i,:).^2) ...
  % ,2));
  %QI = [j j i i]';
  %QJ = [j i j i]';
  %QV = [ ...
  %     Lcot(sub2ind(size(Lcot),i,j)'); ...
  %    -Lcot(sub2ind(size(Lcot),i,j)'); ...
  %    -Lcot(sub2ind(size(Lcot),i,j)'); ...
  %     Lcot(sub2ind(size(Lcot),i,j)')];
  %Q = sparse(QI,QJ,QV,n,n);
  %%Q = -Lcot;
  %Eq = U(:,1)'*Q*U(:,1) + U(:,2)'*Q*U(:,2);
  %LI = [i i];
  %LJ = [ones(3*nt,1) 2*ones(3*nt,1)];
  %LV = [ 2*(repmat(Lcot(sub2ind(size(Lcot),i,j)'),1,2).*(UF(tj,:) - UF(ti,:)) - repmat(Lcot(sub2ind(size(Lcot),k,i)'),1,2).*(UF(ti,:) - UF(tk,:)))];
  %L = sparse(LI(:),LJ(:),LV(:),n,2);
  %El = U(:,1)'*L(:,1)+U(:,2)'*L(:,2);
  %E = E + Eq + El + Ec;
  %[Eq El Ec E]


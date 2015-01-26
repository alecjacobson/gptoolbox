function Q = planarize(V,F,varargin)
  % PLANARIZE Planarize the quads of quad mesh using the local-local technique
  % of "Interactive Planarization and Optimization of 3D Meshes" [Poranne et
  % al. 2012]
  % 
  % Q = planarize(V,F,varargin)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by 4 list of quad indices into V
  %     Optional:
  %       'MaxIter' followed by maximum number of iterations {100}
  %       'mu' followed by "mu" penalty value {0.9}
  %       'r' followed by "r" penalty shrink factor {0.9}
  % Outputs:
  %   Q  #V by dim list of new mesh vertex positions
  %

  % default values
  max_iter = 100;
  mu = 0.9;
  r = 0.9;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','mu','r'}, ...
    {'max_iter','mu','r'});
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

  % number of quads
  m = size(F,1);
  % number of vertices
  n = size(V,1);
  dim = 3;
  T = [F(:,[1 2 3]);F(:,[1 3 4])];
  N = normalizerow(normals(V,T));
  % Initial guess for N
  N = 0.5*(N(1:m,:)+N(m+1:end,:));
  A = doublearea(V,T);
  % A sort of barycentric massmatrix for quad meshes
  M = repdiag(sparse(T(:),T(:),0.25*repmat(A,3,1),n,n),dim);
  E = [];

  % nsp = 3;
  % subplot(1,nsp,1);
  % BC = barycenter(V,F);
  % trisurf(F,V(:,1),V(:,2),V(:,3));
  % hold on;
  % quiver3(BC(:,1),BC(:,2),BC(:,3),N(:,1),N(:,2),N(:,3));
  % hold off;

  % "corner" incidence
  C = sparse(F(:),repmat(1:m,1,4)',1,n,m);

  D = zeros(m,1);

  Q = V;
  iter = 1;
  while true

    % LOCAL
    for ni = 1:m
      Qn = Q(F(ni,:),:);
      Mn = mean(Qn,1);
      Qnc = bsxfun(@minus,Qn,Mn);
      [EV,ED] = eig(Qnc'*Qnc);
      [~,mi] = min(diag(ED));
      N(ni,:) = EV(:,mi);
      D(ni) = -Mn*N(ni,:)';
    end

    % LOCAL
    for v = 1:n
      % faces incident on v
      J = find(C(v,:));
      NJ = N(J,:);
      DJ = D(J);
      AJ = mu*eye(3) + (1-mu)*(NJ'*NJ);
      bJ = mu*V(v,:)' - (1-mu)*NJ'*DJ;
      Q(v,:) = AJ\bJ;
    end

    % GLOBAL
    %% "corner" incidence for each dimension with Nj
    %This is wrong becase D should be scalar...
    %CN = sparse( ...
    %  reshape(bsxfun(@plus,F(:),(0:dim-1)*n),dim*m*4,1), ...
    %  reshape(bsxfun(@plus,repmat((1:m)',4,1),(0:dim-1)*m),dim*m*4,1), ...
    %  reshape(repmat(N,4,1),4*m*dim,1), ...
    %  n*dim,m);
    %% μ ‖P-Q‖² + (1-μ)∑_ij(n_j*q_i-d_j)²
    %% Identity
    %I = repdiag(speye(m,m),dim);
    %H = [          mu*M        (1-mu)*CN; ...
    %          (1-mu)*CN'       (1-mu)*I];
    %l = [-mu*2*M*V(:);zeros(m*dim,1)];
    %QD = min_quad_with_fixed(H,l,[],[]);
    %Q = reshape(QD(1:n*dim),n,dim);
    %D = reshape(QD(n*dim+1:end),m,dim);
    %% (1-μ)∑_ij(n_j*q_i-d_j)²
    %% ‖n_i‖² = 1
    %% 
    %CQ = sparse( ...
    %  reshape(bsxfun(@plus,F(:),(0:dim-1)*n),dim*m*4,1), ...
    %  reshape(bsxfun(@plus,repmat((1:m)',4,1),(0:dim-1)*m),dim*m*4,1), ...
    %  reshape(Q(F(:),:),m*dim*4,1), ...
    %  n*dim,m*dim);
    %%H = (1-mu)*diag(sum(CQ'.^2,2));
    %%l = (1-mu)*-2*diag(sum(CQ',2))*D(:);
    %%% min N'HN - 2N'l
    %%% ‖N‖² = m

    %subplot(1,nsp,2);
    %BC = barycenter(Q,F);
    %trisurf(F,Q(:,1),Q(:,2),Q(:,3));
    %hold on;
    %quiver3(BC(:,1),BC(:,2),BC(:,3),N(:,1),N(:,2),N(:,3));
    %hold off;

    %NT = normalizerow(normals(Q,T));
    %e = sum(sum((NT(1:m,:)-NT(m+(1:m),:)).^2,2));
    %E = [E(:);e];
    %subplot(1,nsp,3);
    %loglog(E);
    iter = iter+1;
    if iter>max_iter
      break;
    end

    mu = r*mu;
  end
end

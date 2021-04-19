function [VV,QQ,SS,J] = catmull_clark(varargin)
  % CATMULL_CLARK Perform iterations of catmull-clark subdivision on a
  % regular, manifold quad mesh (V,Q)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   Q  #Q by 4 list of quad indices into V
  %   Optional:
  %     'Iterations'  followed by number of iterations to conduct {1}
  % Outputs:
  %   VV  #VV by dim list of output mesh vertex positions
  %   QQ  #QQ by 4 list of quad mesh indices into VV
  %   SS  #VV by #V linear subdivision operator so that VV = SS*V
  %   J  #QQ list of indices of birth parents into Q
  % 

  V = varargin{1};
  if nargin>2 && ~ischar(varargin{3})
    poly_input = true;
    PI = varargin{2};
    PC = varargin{3};
    v = 4;
  else
    assert(nargin>=2);
    Q = varargin{2};
    poly_input = false;
    v = 3;
  end

  iters = 1;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Iterations'}, {'iters'});
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

  VV = V;
  SS = speye(size(V,1));
  if poly_input
    QQ = [];
    J = [];
  else
    QQ = Q;
    J = (1:size(Q,1))';
  end

  for i=1:iters
    if i==1 && poly_input
      continue;
    end
    % assume all are quads OR tris
    FF = QQ;
      
    % number of original vertices
    n = size(VV,1);
    % number of original faces
    nf = size(FF,1);

    % extract pure quads
    pure = false(nf,1);
    if size(FF,2)==4
      pure = FF(:,4)~=nan & FF(:,4)>0 & FF(:,4)~=FF(:,3);
      QQ = FF(pure,:);
    else
      QQ = zeros(0,4);
    end
    TT = FF(~pure,1:3);
    JQ = J(pure);
    JT = J(~pure);
    nq = size(QQ,1);
    nt = size(TT,1);
    assert(nf == (nq+nt));
    % face barycenter
    SF = [ ...
      sparse(repmat(1:nq,1,4),QQ(:)',1/4*ones(1,nq*4),nq,n); ...
      sparse(repmat(1:nt,1,3),TT(:)',1/3*ones(1,nt*3),nt,n)];

    % get "quad edges" aka half-edges
    QE = [QQ(:,1:2);QQ(:,2:3);QQ(:,3:4);QQ(:,[4 1])];
    % get "triangle edges"
    TE = [TT(:,1:2);TT(:,2:3);TT(:,[3 1])];
    % all edges
    FE = [QE;TE];
    % get unique edge list, QE2E tells original face edges (with duplicates,
    % EE) where to find matching unique edge (E)
    [E,I,FE2E] = unique(sort(FE,2),'rows');
    FEF = [repmat(1:nq,1,4) nq+repmat(1:nt,1,3)]';
    % number of original edges
    ne = size(E,1);
    E2F = zeros(ne,2);
    flip = FE(:,1)~=E(FE2E,1);
    E2F(sub2ind(size(E2F),FE2E,flip+1)) = FEF;
    Eout = accumarray(FE2E,1,[size(E,1) 1])==1;
    in = find(~Eout);
    out = find(Eout);
    % Average of original edge points and edge-flap face barycenters
    SE = ...
      (sparse( ...
        repmat(in,4,1), ...
        [E(in,1);E(in,2);n+[E2F(in,1);E2F(in,2)]], ...
        1/4*ones(numel(in)*4,1),ne,n+nf) + ...
      sparse( ...
        repmat(out,2,1), ...
        [E(out,1);E(out,2)], ...
        repmat(1/2,numel(out)*2,1),ne,n+nf)) ...
      *[speye(n);SF];


    V2F = [ ...
      sparse(QQ,repmat(1:nq,4,1)',1,n,nq) ...
      sparse(TT,repmat(1:nt,3,1)',1,n,nt)];
    % number of faces incident on each vertex
    val = sum(V2F,2);
    % http://www.alecjacobson.com/weblog/?p=3235
    pou_rows =  @(A) spdiags (1./sum (A,2), 0, size(A,1), size(A,1)) * A ;
    % normalize to take average
    %V2F = bsxfun(@rdivide,V2F,sum(V2F,2));
    V2F = pou_rows(V2F);
    % compute midpoints of original vertices
    ME = sparse(E(:,[1 2 1 2]),E(:,[1 2 2 1]),1,n,n);
    %ME = bsxfun(@rdivide,ME,sum(ME,2));
    ME = pou_rows(ME);
    % Weighting 
    W = bsxfun(@rdivide,[(val-3) ones(n,1) 2*ones(n,1)],val);
    Vout = false(n,1);
    Vout(E(Eout,:)) = true;
    in = find(~Vout);
    out = find(Vout);
    % 1/8-------3/4-------1/8
    %    1/2--1/2  1/2--1/2
    %       1/4 1/2 1/4
    SV = ( ...
      sparse(repmat(in,1,3),in+(0:2)*n, W(in,:),n,3*n) + ...
      sparse(E(Eout,[1 2 1 2]),E(Eout,[1 2 2 1]),0.125,n,3*n) + ...
      sparse(out,out,0.5,n,3*n)) ...
      *[speye(n);V2F*SF;ME];

    S = [SV;SF;SE];
    VV = S*VV;

    %clf;
    %hold on;
    %tsurf(Q,V,'FaceAlpha',0.5);
    %scatter3(VV(:,1),VV(:,2),VV(:,3),'LineWidth',3);
    %hold off;
    %error

    SS = S*SS;
    QE2E = reshape(FE2E(1:nq*4),size(QQ));
    TE2E = reshape(FE2E(nq*4+(1:nt*3)),size(TT));
    QQ = [ ...
      QQ(:,1) n+nf+QE2E(:,1) n+(1:nq)' n+nf+QE2E(:,4); ...
      QQ(:,2) n+nf+QE2E(:,2) n+(1:nq)' n+nf+QE2E(:,1); ...
      QQ(:,3) n+nf+QE2E(:,3) n+(1:nq)' n+nf+QE2E(:,2); ...
      QQ(:,4) n+nf+QE2E(:,4) n+(1:nq)' n+nf+QE2E(:,3); ...
      ];
    TT = [ ...
      TT(:,1) n+nf+TE2E(:,1) n+nq+(1:nt)' n+nf+TE2E(:,3); ...
      TT(:,2) n+nf+TE2E(:,2) n+nq+(1:nt)' n+nf+TE2E(:,1); ...
      TT(:,3) n+nf+TE2E(:,3) n+nq+(1:nt)' n+nf+TE2E(:,2)];
    QQ = [QQ;TT];
    J = [JQ;JQ;JQ;JQ;JT;JT;JT];

  end
      
    
end
